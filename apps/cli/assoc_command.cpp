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

#include "data/loader/bim_loader.h"
#include "estimator/bayes/indicator.h"
#include "gelex/types/assoc_input.h"
#include "utils/formatter.h"
#include "utils/utils.h"

auto assoc_command(argparse::ArgumentParser& cmd) -> void
{
    cmd.add_description(
        "Perform genome-wide association study using mixed linear model");

    // ================================================================
    // IO
    // ================================================================
    cmd.add_group("Data Files");
    cmd.add_argument("-p", "--pheno")
        .help("Phenotype file (TSV format: FID, IID, trait1, ...)")
        .metavar("<PHENOTYPE>")
        .required();
    cmd.add_argument("--pheno-col")
        .help("Phenotype column index (0-based)")
        .default_value(2)
        .scan<'i', int>();
    cmd.add_argument("-b", "--bfile")
        .help("PLINK binary file prefix (.bed/.bim/.fam)")
        .metavar("<BFILE>")
        .required();
    cmd.add_argument("--grm")
        .help("GRM file prefix(es). Can specify multiple GRMs.")
        .metavar("<GRM>")
        .nargs(argparse::nargs_pattern::at_least_one)
        .required();
    cmd.add_argument("--qcovar")
        .default_value("")
        .help("Quantitative covariates (TSV: FID, IID, covar1, ...)");
    cmd.add_argument("--dcovar")
        .default_value("")
        .help("Discrete covariates (TSV: FID, IID, factor1, ...)");
    cmd.add_argument("-o", "--out")
        .help("Output file prefix")
        .metavar("<OUT>")
        .required()
        .default_value("gelex");
    // ================================================================
    // REML Configuration
    // ================================================================
    cmd.add_group("REML Options");
    cmd.add_argument("--max-iter")
        .help("Max iteration in REML process")
        .default_value(100)
        .scan<'i', int>();
    cmd.add_argument("--tol")
        .help("tolerance for convergence in REML process")
        .default_value(1e-6)
        .scan<'g', double>();

    // ================================================================
    // Data Processing
    // ================================================================
    cmd.add_group("Processing Options");

    cmd.add_argument("--chunk-size")
        .help("SNPs per chunk for association testing")
        .default_value(1000)
        .scan<'i', int>();
    cmd.add_argument("--iid-only")
        .help("Use only IID for sample matching (ignore FID)")
        .flag();
    cmd.add_argument("--loco")
        .help("Enable Leave-One-Chromosome-Out (LOCO) mode")
        .flag();

    // ================================================================
    // Model Configuration
    // ================================================================
    cmd.add_group("Model Configuration");
    cmd.add_argument("--model")
        .help("Association model: a (additive), d (dominance), ad (both)")
        .default_value("a")
        .metavar("<MODEL>")
        .choices("a", "d", "ad");
    // ================================================================
    // Performance
    // ================================================================
    cmd.add_group("Performance");
    const int default_threads
        = std::max(1U, std::thread::hardware_concurrency() / 2);
    cmd.add_argument("--threads")
        .help("Number of CPU threads to use")
        .default_value(default_threads)
        .scan<'i', int>();
}

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
    auto snp_effects
        = gelex::detail::BimLoader(bed_path.replace_extension(".bim"))
              .take_info();

    const Eigen::Index n_snps = bed_pipe.num_snps();
    const auto n_samples
        = static_cast<Eigen::Index>(sample_manager->num_common_samples());
    const int chunk_size = cmd.get<int>("--chunk-size");

    // Group SNPs by chromosome
    struct ChrRange
    {
        std::string chrom;
        size_t start_idx;
        size_t end_idx;
    };
    std::vector<ChrRange> chr_ranges;
    if (n_snps > 0)
    {
        std::string current_chr = snp_effects[0].chrom;
        size_t start = 0;
        for (size_t i = 1; i < static_cast<size_t>(n_snps); ++i)
        {
            if (snp_effects[i].chrom != current_chr)
            {
                chr_ranges.push_back({current_chr, start, i});
                current_chr = snp_effects[i].chrom;
                start = i;
            }
        }
        chr_ranges.push_back({current_chr, start, static_cast<size_t>(n_snps)});
    }

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

    // barkeep progress bar
    size_t progress_counter{0};
    auto pbar = gelex::detail::create_association_progress_bar(
        progress_counter, static_cast<size_t>(n_snps));
    pbar.pbar->show();

    auto run_assoc_for_range = [&](const ChrRange& range,
                                   Eigen::MatrixXd&& v_inv,
                                   Eigen::VectorXd&& v_inv_residual)
    {
        gelex::AssocInput input(
            chunk_size, std::move(v_inv), std::move(v_inv_residual));
        gelex::AssocOutput output(chunk_size);
        Eigen::VectorXd freqs(chunk_size);

        size_t range_len = range.end_idx - range.start_idx;
        size_t n_chunks = (range_len + chunk_size - 1) / chunk_size;

        for (size_t chunk_idx = 0; chunk_idx < n_chunks; ++chunk_idx)
        {
            size_t start = range.start_idx + (chunk_idx * chunk_size);
            size_t end = std::min(start + chunk_size, range.end_idx);
            Eigen::Index current_chunk_size = end - start;

            input.Z.resize(n_samples, current_chunk_size);
            freqs.resize(current_chunk_size);

            bed_pipe.load_chunk(
                input.Z,
                static_cast<Eigen::Index>(start),
                static_cast<Eigen::Index>(end));
            encoder(input.Z, additive, &freqs);
            gelex::gwas::wald_test(input, output);

            for (size_t i = 0; i < static_cast<size_t>(current_chunk_size); ++i)
            {
                writer.write_result(
                    snp_effects[start + i],
                    {.freq = freqs(static_cast<Eigen::Index>(i)),
                     .beta = output.beta(static_cast<Eigen::Index>(i)),
                     .se = output.se(static_cast<Eigen::Index>(i)),
                     .p_value = output.p_value(static_cast<Eigen::Index>(i))});
            }

            progress_counter += static_cast<size_t>(current_chunk_size);
            double finish_percent = static_cast<double>(progress_counter)
                                    / static_cast<double>(n_snps) * 100;

            pbar.status->message(
                fmt::format(
                    "{:.1f}% ({}/{}) | {}",
                    finish_percent,
                    gelex::HumanReadable(progress_counter),
                    gelex::HumanReadable(n_snps),
                    eta_calculator.get_eta(progress_counter)));
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
        gelex::Estimator estimator(
            cmd.get<int>("--max-iter"), cmd.get<double>("--tol"));

        if (model.genetic().size() != loco_loaders.size())
        {
            throw gelex::InvalidInputException(
                "Number of genetic components in model does not match number "
                "of GRMs provided.");
        }

        for (const auto& range : chr_ranges)
        {
            for (size_t i = 0; i < loco_loaders.size(); ++i)
            {
                // Load LOCO GRM for this chromosome for each loader
                std::filesystem::path chr_grm_prefix = grm_paths[i];
                chr_grm_prefix += ".chr" + range.chrom;

                loco_loaders[i].load_loco_grm(
                    chr_grm_prefix,
                    sample_manager->common_id_map(),
                    model.genetic()[i].K);
            }

            // Re-fit model
            Eigen::MatrixXd v_inv
                = estimator.fit(model, state, true, false);  // verbose=false
            Eigen::VectorXd v_inv_residual
                = v_inv
                  * (model.phenotype() - model.fixed().X * state.fixed().coeff);

            run_assoc_for_range(
                range, std::move(v_inv), std::move(v_inv_residual));
        }
    }
    else
    {
        // Standard mode: one REML fit for all chromosomes
        gelex::FreqModel model(data_pipe);
        gelex::FreqState state(model);
        gelex::Estimator estimator(
            cmd.get<int>("--max-iter"), cmd.get<double>("--tol"));

        Eigen::MatrixXd v_inv = estimator.fit(model, state, true, true);
        Eigen::VectorXd v_inv_residual
            = v_inv
              * (model.phenotype() - model.fixed().X * state.fixed().coeff);

        for (const auto& range : chr_ranges)
        {
            run_assoc_for_range(
                range, std::move(v_inv), std::move(v_inv_residual));
        }
    }

    pbar.pbar->done();
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
