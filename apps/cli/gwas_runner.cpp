#include "cli/gwas_runner.h"

#include <Eigen/Core>
#include <functional>
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
      chr_groups_(build_chr_groups(config_.loco, snp_effects_)),
      assoc_input_(
          static_cast<Eigen::Index>(
              data_pipe_.sample_manager()->num_common_samples()),
          static_cast<Eigen::Index>(config_.chunk_size)),
      assoc_output_(static_cast<Eigen::Index>(config_.chunk_size)),
      freqs_(static_cast<Eigen::Index>(config_.chunk_size))
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

auto GwasRunner::update_assoc_input(
    const FreqModel& model,
    const FreqState& state,
    Eigen::MatrixXd&& v_inv) -> void
{
    assoc_input_.V_inv = std::move(v_inv);
    assoc_input_.V_inv_y
        = assoc_input_.V_inv
          * (model.phenotype() - model.fixed().X * state.fixed().coeff);
}

auto GwasRunner::run_normal() -> void
{
    FreqModel model(data_pipe_);
    FreqState state(model);
    Estimator estimator(config_.max_iter, config_.tol);

    update_assoc_input(model, state, estimator.fit(model, state, true, true));

    std::atomic<size_t> progress_counter{0};
    auto pbar = create_progress_bar(progress_counter, snp_effects_.size());
    pbar.display->show();

    auto progress_callback = [&](size_t current, size_t total, size_t offset)
    {
        pbar.status->message(
            fmt::format(
                "{:.1f}% ({}/{}) | {}",
                static_cast<double>(current) / static_cast<double>(total)
                    * 100.0,
                HumanReadable(current),
                HumanReadable(total),
                eta_calculator_.get_eta(current + offset)));
    };

    for (const auto& group : chr_groups_)
    {
        scan_chromosome(
            group, progress_counter, snp_effects_.size(), 0, progress_callback);
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

    const auto id_map = data_pipe_.sample_manager()->common_id_map();
    std::vector<gelex::LocoGRMLoader> loco_loaders;
    loco_loaders.reserve(config_.grm_paths.size());
    for (const auto& path : config_.grm_paths)
    {
        loco_loaders.emplace_back(path, id_map);
    }

    size_t total_processed = 0;

    for (const auto& group : chr_groups_)
    {
        for (size_t i = 0; i < loco_loaders.size(); ++i)
        {
            const auto chr_grm_prefix
                = config_.grm_paths[i].string() + ".chr" + group.name;
            loco_loaders[i].load_loco_grm(
                chr_grm_prefix, id_map, model.genetic()[i].K);
        }

        auto loco_logger
            = std::make_unique<gelex::detail::LocoRemlLogger>(group.name);
        Estimator estimator(
            config_.max_iter, config_.tol, std::move(loco_logger));

        update_assoc_input(
            model, state, estimator.fit(model, state, true, true));

        std::atomic<size_t> chr_counter{0};
        auto chr_pbar = create_progress_bar(chr_counter, group.total_snps);
        chr_pbar.display->show();

        auto progress_callback
            = [&](size_t current, size_t total, size_t offset)
        {
            chr_pbar.status->message(
                fmt::format(
                    "{:.1f}% ({}/{}) | {}",
                    static_cast<double>(current) / static_cast<double>(total)
                        * 100.0,
                    HumanReadable(current),
                    HumanReadable(total),
                    eta_calculator_.get_eta(offset + current)));
        };

        scan_chromosome(
            group,
            chr_counter,
            group.total_snps,
            total_processed,
            progress_callback);

        total_processed += group.total_snps;
        chr_pbar.display->done();
    }
}

auto GwasRunner::scan_chromosome(
    const ChrGroup& group,
    std::atomic<size_t>& progress_counter,
    size_t total_snps_to_report,
    size_t total_processed_before,
    const std::function<void(size_t, size_t, size_t)>& progress_callback)
    -> void
{
    const auto n_samples = assoc_input_.V_inv.rows();
    const auto chunk_size = static_cast<size_t>(config_.chunk_size);

    for (const auto& [range_start, range_end] : group.ranges)
    {
        for (auto start = range_start; start < range_end; start += chunk_size)
        {
            const auto end
                = std::min(start + chunk_size, static_cast<size_t>(range_end));
            const auto current_chunk_size
                = static_cast<Eigen::Index>(end - start);

            assoc_input_.resize(n_samples, current_chunk_size);
            assoc_output_.resize(current_chunk_size);
            freqs_.resize(current_chunk_size);

            bed_pipe_.load_chunk(
                assoc_input_.Z,
                static_cast<Eigen::Index>(start),
                static_cast<Eigen::Index>(end));

            encoder_(assoc_input_.Z, config_.additive, &freqs_);
            gwas::wald_test(assoc_input_, assoc_output_);

            for (Eigen::Index i = 0; i < current_chunk_size; ++i)
            {
                writer_.write_result(
                    snp_effects_[start + static_cast<size_t>(i)],
                    {.freq = freqs_(i),
                     .beta = assoc_output_.beta(i),
                     .se = assoc_output_.se(i),
                     .p_value = assoc_output_.p_value(i)});
            }

            const auto current_progress
                = progress_counter.fetch_add(
                      static_cast<size_t>(current_chunk_size),
                      std::memory_order_relaxed)
                  + static_cast<size_t>(current_chunk_size);

            if (progress_callback)
            {
                progress_callback(
                    current_progress,
                    total_snps_to_report,
                    total_processed_before);
            }
        }
    }
}

}  // namespace gelex::cli
