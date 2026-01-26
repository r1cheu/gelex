#include "assoc_command.h"

#include <filesystem>
#include <memory>
#include <ranges>
#include <thread>

#include <fmt/format.h>
#include <Eigen/Core>

#include "cli_helper.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/data_pipe.h"
#include "gelex/data/grm_code_policy.h"
#include "gelex/data/loco_grm_loader.h"
#include "gelex/estimator/freq/estimator.h"
#include "gelex/estimator/freq/reml.h"
#include "gelex/gwas/association_test.h"
#include "gelex/gwas/gwas_writer.h"
#include "gelex/logger.h"
#include "gelex/model/freq/model.h"
#include "logger/loco_reml_logger.h"

#include "data/loader/bim_loader.h"
#include "gelex/types/assoc_input.h"
#include "utils/formatter.h"
#include "utils/utils.h"

auto assoc_execute(argparse::ArgumentParser& cmd) -> int
{
    auto logger = gelex::logging::get();
    std::string out_prefix = cmd.get("--out");

    // Parse model and test type
    // auto assoc_mode = gelex::gwas::parse_assoc_mode(cmd.get("--model"));

    gelex::cli::setup_parallelization(cmd.get<int>("--threads"));

    auto grm_paths = std::ranges::to<std::vector<std::filesystem::path>>(
        cmd.get<std::vector<std::string>>("--grm"));

    auto bed_path = gelex::BedPipe::format_bed_path(cmd.get("bfile"));

    // ================================================================
    // Data Loading
    // ================================================================
    gelex::DataPipe::Config config{
        .phenotype_path = cmd.get("pheno"),
        .phenotype_column = cmd.get<int>("--pheno-col"),
        .bed_path = bed_path,
        .use_dominance_effect = false,
        .use_mmap = false,
        .chunk_size = cmd.get<int>("--chunk-size"),
        .qcovar_path = cmd.get("--qcovar"),
        .dcovar_path = cmd.get("--dcovar"),
        .iid_only = cmd.get<bool>("--iid-only"),
        .output_prefix = cmd.get("--out"),
        .grm_paths = grm_paths};

    gelex::cli::print_assoc_header(cmd.get<int>("--threads"));

    auto data_pipe = gelex::load_data_for_reml(config);
    auto sample_manager = data_pipe.sample_manager();

    bool loco = cmd.get<bool>("--loco");

    gelex::BedPipe bed_pipe(bed_path, sample_manager);
    auto bim_path = bed_path;
    bim_path.replace_extension(".bim");
    auto snp_effects = gelex::detail::BimLoader(bim_path).take_info();

    const Eigen::Index n_snps = bed_pipe.num_snps();
    const auto n_samples
        = static_cast<Eigen::Index>(sample_manager->num_common_samples());
    const int chunk_size = cmd.get<int>("--chunk-size");

    auto chr_groups = gelex::cli::build_chr_groups(true, snp_effects);

    logger->info("");
    logger->info(gelex::section("Running Association Tests..."));
    logger->info(gelex::task("SNPs to test : {}", n_snps));
    logger->info(gelex::task("Chunk size   : {}", chunk_size));
    if (loco)
    {
        logger->info(gelex::task("Mode         : LOCO"));
    }
    logger->info("");

    gelex::gwas::GwasWriter writer(out_prefix);
    writer.write_header();

    gelex::grm::Yang encoder;
    bool additive = cmd.get("--model") == "a";
    gelex::SmoothEtaCalculator eta_calculator(n_snps);

    std::atomic<size_t> progress_counter{0};
    auto pbar = gelex::cli::create_progress_bar(
        progress_counter, static_cast<size_t>(n_snps));

    auto process_group
        = [&](gelex::AssocInput& input, const gelex::cli::ChrGroup& group)
    {
        gelex::AssocOutput output(chunk_size);
        Eigen::VectorXd freqs(chunk_size);

        for (const auto& [range_start, range_end] : group.ranges)
        {
            auto range_len = static_cast<size_t>(range_end - range_start);
            auto n_chunks = (range_len + chunk_size - 1) / chunk_size;

            for (size_t chunk_idx = 0; chunk_idx < n_chunks; ++chunk_idx)
            {
                auto start = range_start + (chunk_idx * chunk_size);
                auto end = std::min(
                    start + chunk_size, static_cast<size_t>(range_end));
                auto current_chunk_size
                    = static_cast<Eigen::Index>(end - start);

                input.Z.resize(n_samples, current_chunk_size);
                freqs.resize(current_chunk_size);

                bed_pipe.load_chunk(
                    input.Z,
                    static_cast<Eigen::Index>(start),
                    static_cast<Eigen::Index>(end));
                encoder(input.Z, additive, &freqs);
                gelex::gwas::wald_test(input, output);

                for (Eigen::Index i = 0; i < current_chunk_size; ++i)
                {
                    writer.write_result(
                        snp_effects[start + static_cast<size_t>(i)],
                        {.freq = freqs(i),
                         .beta = output.beta(i),
                         .se = output.se(i),
                         .p_value = output.p_value(i)});
                }

                auto current_progress = progress_counter.fetch_add(
                    static_cast<size_t>(current_chunk_size),
                    std::memory_order_relaxed);
                current_progress += static_cast<size_t>(current_chunk_size);

                pbar.status->message(
                    fmt::format(
                        "{:.1f}% ({}/{}) | {}",
                        static_cast<double>(current_progress) / n_snps * 100,
                        gelex::HumanReadable(current_progress),
                        gelex::HumanReadable(n_snps),
                        eta_calculator.get_eta(current_progress)));
            }
        }
    };

    if (loco)
    {
        std::vector<gelex::LocoGRMLoader> loco_loaders;
        loco_loaders.reserve(grm_paths.size());
        for (const auto& path : grm_paths)
        {
            loco_loaders.emplace_back(path, sample_manager->common_id_map());
        }

        gelex::FreqModel model(data_pipe);
        gelex::FreqState state(model);

        if (model.genetic().size() != loco_loaders.size())
        {
            throw gelex::InvalidInputException(
                "Number of genetic components in model does not match number "
                "of GRMs provided.");
        }

        for (const auto& group : chr_groups)
        {
            for (size_t i = 0; i < loco_loaders.size(); ++i)
            {
                auto chr_grm_prefix = grm_paths[i];
                chr_grm_prefix += ".chr" + group.name;

                loco_loaders[i].load_loco_grm(
                    chr_grm_prefix,
                    sample_manager->common_id_map(),
                    model.genetic()[i].K);
            }

            auto loco_logger
                = std::make_unique<gelex::detail::LocoRemlLogger>(group.name);
            gelex::Estimator estimator(
                cmd.get<int>("--max-iter"),
                cmd.get<double>("--tol"),
                std::move(loco_logger));

            auto v_inv = estimator.fit(model, state, true, true);
            auto v_inv_residual
                = v_inv
                  * (model.phenotype() - model.fixed().X * state.fixed().coeff);

            // Create per-chromosome progress bar
            std::atomic<size_t> chr_counter{0};
            auto chr_pbar = gelex::cli::create_progress_bar(
                chr_counter, static_cast<size_t>(group.total_snps));
            chr_pbar.display->show();

            gelex::AssocInput input(
                chunk_size, std::move(v_inv), std::move(v_inv_residual));

            // Inline scan for this chromosome with local progress
            gelex::AssocOutput output(chunk_size);
            Eigen::VectorXd freqs(chunk_size);

            for (const auto& [range_start, range_end] : group.ranges)
            {
                auto range_len = static_cast<size_t>(range_end - range_start);
                auto n_chunks = (range_len + chunk_size - 1) / chunk_size;

                for (size_t chunk_idx = 0; chunk_idx < n_chunks; ++chunk_idx)
                {
                    auto start = range_start + (chunk_idx * chunk_size);
                    auto end = std::min(
                        start + chunk_size, static_cast<size_t>(range_end));
                    auto current_chunk_size
                        = static_cast<Eigen::Index>(end - start);

                    input.Z.resize(n_samples, current_chunk_size);
                    freqs.resize(current_chunk_size);

                    bed_pipe.load_chunk(
                        input.Z,
                        static_cast<Eigen::Index>(start),
                        static_cast<Eigen::Index>(end));
                    encoder(input.Z, additive, &freqs);
                    gelex::gwas::wald_test(input, output);

                    for (Eigen::Index j = 0; j < current_chunk_size; ++j)
                    {
                        writer.write_result(
                            snp_effects[start + static_cast<size_t>(j)],
                            {.freq = freqs(j),
                             .beta = output.beta(j),
                             .se = output.se(j),
                             .p_value = output.p_value(j)});
                    }

                    auto chr_progress = chr_counter.fetch_add(
                        static_cast<size_t>(current_chunk_size),
                        std::memory_order_relaxed);
                    chr_progress += static_cast<size_t>(current_chunk_size);

                    progress_counter.fetch_add(
                        static_cast<size_t>(current_chunk_size),
                        std::memory_order_relaxed);

                    chr_pbar.status->message(
                        fmt::format(
                            "{:.1f}% ({}/{}) | {}",
                            static_cast<double>(chr_progress) / group.total_snps
                                * 100,
                            gelex::HumanReadable(chr_progress),
                            gelex::HumanReadable(group.total_snps),
                            eta_calculator.get_eta(progress_counter.load())));
                }
            }
            chr_pbar.display->done();
        }
    }
    else
    {
        gelex::FreqModel model(data_pipe);
        gelex::FreqState state(model);
        gelex::Estimator estimator(
            cmd.get<int>("--max-iter"), cmd.get<double>("--tol"));

        auto v_inv = estimator.fit(model, state, true, true);
        auto v_inv_residual
            = v_inv
              * (model.phenotype() - model.fixed().X * state.fixed().coeff);

        pbar.display->show();

        gelex::AssocInput input(
            chunk_size, std::move(v_inv), std::move(v_inv_residual));
        for (const auto& group : chr_groups)
        {
            process_group(input, group);
        }

        pbar.display->done();
    }
    writer.finalize();

    logger->info("");
    logger->info(
        gelex::success(
            "Scan complete! Time elapsed: {}",
            eta_calculator.total_time_consumed()));
    logger->info(gelex::success("Results saved to : {}.gwas.tsv", out_prefix));
    logger->info(
        fmt::format(
            fmt::emphasis::bold | fmt::fg(fmt::color::light_cyan),
            "───────────────────────────────────"
            "───────────────────────────────────"));

    return 0;
}
