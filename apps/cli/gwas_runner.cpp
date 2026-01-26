#include "cli/gwas_runner.h"

#include <Eigen/Core>
#include <memory>

#include <fmt/format.h>

#include "gelex/data/grm_code_policy.h"
#include "gelex/estimator/freq/estimator.h"
#include "gelex/gwas/association_test.h"
#include "gelex/logger.h"
#include "gelex/model/freq/model.h"
#include "logger/loco_reml_logger.h"
#include "utils/formatter.h"

namespace gelex::cli
{

GwasRunner::GwasRunner(
    GwasRunner::Config config,
    DataPipe data_pipe,
    BedPipe bed_pipe,
    SnpEffects snp_effects)
    : config_(std::move(config)),
      data_pipe_(std::move(data_pipe)),
      bed_pipe_(std::move(bed_pipe)),
      writer_(config_.out_prefix),
      snp_effects_(std::move(snp_effects)),
      eta_calculator_(static_cast<Eigen::Index>(snp_effects_.size())),
      chr_groups_(build_chr_groups(config_.loco, snp_effects_))
{
}

auto GwasRunner::run() -> void
{
    print_summary();
    writer_.write_header();

    if (config_.loco)
    {
        run_loco();
    }
    else
    {
        run_normal();
    }

    writer_.finalize();

    auto logger = gelex::logging::get();
    logger->info("");
    logger->info(
        gelex::success(
            "Scan complete! Time elapsed: {}",
            eta_calculator_.total_time_consumed()));
    logger->info(
        gelex::success("Results saved to : {}.gwas.tsv", config_.out_prefix));
    logger->info(
        fmt::format(
            fmt::emphasis::bold | fmt::fg(fmt::color::light_cyan),
            "───────────────────────────────────"
            "───────────────────────────────────"));
}

auto GwasRunner::print_summary() const -> void
{
    auto logger = gelex::logging::get();
    logger->info("");
    logger->info(gelex::section("Running Association Tests..."));
    logger->info(gelex::task("SNPs to test : {}", snp_effects_.size()));
    logger->info(gelex::task("Chunk size   : {}", config_.chunk_size));

    if (config_.loco)
    {
        logger->info(gelex::task("Mode         : LOCO"));
    }
    logger->info("");
}

auto GwasRunner::run_normal() -> void
{
    FreqModel model(data_pipe_);
    FreqState state(model);
    Estimator estimator(config_.max_iter, config_.tol);

    auto v_inv = estimator.fit(model, state, true, true);
    Eigen::VectorXd v_inv_residual
        = v_inv * (model.phenotype() - model.fixed().X * state.fixed().coeff);

    std::atomic<size_t> progress_counter{0};
    auto pbar = create_progress_bar(progress_counter, snp_effects_.size());
    pbar.display->show();

    AssocInput input(
        config_.chunk_size, std::move(v_inv), std::move(v_inv_residual));

    auto progress_callback = [&](size_t current, size_t total)
    {
        pbar.status->message(
            fmt::format(
                "{:.1f}% ({}/{}) | {}",
                static_cast<double>(current) / static_cast<double>(total)
                    * 100.0,
                HumanReadable(current),
                HumanReadable(total),
                eta_calculator_.get_eta(current)));
    };

    for (const auto& group : chr_groups_)
    {
        scan_chromosome(
            input,
            group,
            progress_counter,
            snp_effects_.size(),
            progress_callback);
    }

    pbar.display->done();
}

auto GwasRunner::run_loco() -> void
{
    FreqModel model(data_pipe_);
    FreqState state(model);

    if (model.genetic().size() != config_.grm_paths.size())
    {
        throw InvalidInputException(
            "Number of genetic components in model does not match number "
            "of GRMs provided.");
    }

    std::vector<gelex::LocoGRMLoader> loco_loaders;
    loco_loaders.reserve(config_.grm_paths.size());
    auto id_map = data_pipe_.sample_manager()->common_id_map();
    for (const auto& path : config_.grm_paths)
    {
        loco_loaders.emplace_back(path, id_map);
    }

    std::atomic<size_t> total_progress_counter{0};

    for (const auto& group : chr_groups_)
    {
        for (size_t i = 0; i < loco_loaders.size(); ++i)
        {
            auto chr_grm_prefix = config_.grm_paths[i];
            chr_grm_prefix += ".chr" + group.name;

            loco_loaders[i].load_loco_grm(
                chr_grm_prefix,
                data_pipe_.sample_manager()->common_id_map(),
                model.genetic()[i].K);
        }

        auto loco_logger
            = std::make_unique<gelex::detail::LocoRemlLogger>(group.name);
        Estimator estimator(
            config_.max_iter, config_.tol, std::move(loco_logger));

        auto v_inv = estimator.fit(model, state, true, true);
        Eigen::VectorXd v_inv_residual
            = v_inv
              * (model.phenotype() - model.fixed().X * state.fixed().coeff);

        std::atomic<size_t> chr_counter{0};
        auto chr_pbar = create_progress_bar(chr_counter, group.total_snps);
        chr_pbar.display->show();

        AssocInput input(
            config_.chunk_size, std::move(v_inv), std::move(v_inv_residual));

        auto progress_callback = [&](size_t current, size_t total)
        {
            chr_pbar.status->message(
                fmt::format(
                    "{:.1f}% ({}/{}) | {}",
                    static_cast<double>(current) / static_cast<double>(total)
                        * 100.0,
                    HumanReadable(current),
                    HumanReadable(total),
                    eta_calculator_.get_eta(total_progress_counter.load())));
        };

        scan_chromosome(
            input, group, chr_counter, group.total_snps, progress_callback);

        total_progress_counter.fetch_add(group.total_snps);
        chr_pbar.display->done();
    }
}

auto GwasRunner::scan_chromosome(
    AssocInput& input,
    const ChrGroup& group,
    std::atomic<size_t>& progress_counter,
    size_t total_snps,
    const std::function<void(size_t, size_t)>& progress_callback) -> void
{
    const Eigen::Index n_samples = input.V_inv.rows();
    AssocOutput output(config_.chunk_size);
    Eigen::VectorXd freqs(config_.chunk_size);

    for (const auto& [range_start, range_end] : group.ranges)
    {
        auto range_len = static_cast<size_t>(range_end - range_start);
        auto n_chunks
            = (range_len + config_.chunk_size - 1) / config_.chunk_size;

        for (size_t chunk_idx = 0; chunk_idx < n_chunks; ++chunk_idx)
        {
            auto start = range_start + (chunk_idx * config_.chunk_size);
            auto end = std::min(
                start + config_.chunk_size, static_cast<size_t>(range_end));
            auto current_chunk_size = static_cast<Eigen::Index>(end - start);

            input.Z.resize(n_samples, current_chunk_size);
            freqs.resize(current_chunk_size);

            bed_pipe_.load_chunk(
                input.Z,
                static_cast<Eigen::Index>(start),
                static_cast<Eigen::Index>(end));
            encoder_(input.Z, config_.additive, &freqs);
            gwas::wald_test(input, output);

            for (Eigen::Index i = 0; i < current_chunk_size; ++i)
            {
                writer_.write_result(
                    snp_effects_[start + static_cast<size_t>(i)],
                    {.freq = freqs(i),
                     .beta = output.beta(i),
                     .se = output.se(i),
                     .p_value = output.p_value(i)});
            }

            auto current_progress = progress_counter.fetch_add(
                static_cast<size_t>(current_chunk_size),
                std::memory_order_relaxed);
            current_progress += static_cast<size_t>(current_chunk_size);

            if (progress_callback)
            {
                progress_callback(current_progress, total_snps);
            }
        }
    }
}

}  // namespace gelex::cli
