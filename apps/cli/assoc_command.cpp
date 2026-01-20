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
#include "gelex/estimator/freq/reml.h"
#include "gelex/gwas/association_test.h"
#include "gelex/gwas/gwas_writer.h"
#include "gelex/logger.h"

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
    auto [sample_manager, v_inv, v_inv_residual] = gelex::reml(
        config,
        cmd.get<int>("--max-iter"),
        cmd.get<double>("--tol"),
        true,
        true);

    gelex::BedPipe bed_pipe(bed_path, sample_manager);
    auto snp_effects
        = gelex::detail::BimLoader(bed_path.replace_extension(".bim"))
              .take_info();

    const Eigen::Index n_snps = bed_pipe.num_snps();
    const auto n_samples
        = static_cast<Eigen::Index>(sample_manager->num_common_samples());
    const int chunk_size = cmd.get<int>("--chunk-size");

    logger->info("");
    logger->info(gelex::section("Running Association Tests..."));
    logger->info(gelex::task("SNPs to test : {}", n_snps));
    logger->info(gelex::task("Chunk size   : {}", chunk_size));
    logger->info("");

    gelex::AssocInput input(
        chunk_size, std::move(v_inv), std::move(v_inv_residual));
    gelex::AssocOutput output(chunk_size);

    // ================================================================
    // Association testing
    // ================================================================

    gelex::gwas::GwasWriter writer(out_prefix);
    writer.write_header();
    size_t n_tested = 0;
    auto n_chunks = static_cast<size_t>((n_snps + chunk_size - 1) / chunk_size);

    // barkeep progress bar
    size_t progress_counter{0};
    auto pbar = gelex::detail::create_association_progress_bar(
        progress_counter, static_cast<size_t>(n_snps));
    pbar.pbar->show();

    gelex::grm::Yang encoder;
    bool additive = cmd.get("--model") == "a";

    gelex::SmoothEtaCalculator eta_calculator(n_snps);

    for (size_t chunk_idx = 0; chunk_idx < n_chunks; ++chunk_idx)
    {
        Eigen::Index start = static_cast<Eigen::Index>(chunk_idx)
                             * static_cast<Eigen::Index>(chunk_size);
        Eigen::Index end
            = std::min(start + static_cast<Eigen::Index>(chunk_size), n_snps);
        input.Z.resize(n_samples, (end - start));
        bed_pipe.load_chunk(input.Z, start, end);
        encoder(input.Z, additive);
        gelex::gwas::wald_test(input, output);

        for (size_t i = n_tested;
             i < n_tested + static_cast<size_t>(end - start);
             ++i)
        {
            writer.write_result(
                snp_effects[i],
                {.beta = output.beta(static_cast<Eigen::Index>(i - n_tested)),
                 .se = output.se(static_cast<Eigen::Index>(i - n_tested)),
                 .p_value
                 = output.p_value(static_cast<Eigen::Index>(i - n_tested))});
        }

        n_tested += static_cast<size_t>(end - start);
        progress_counter += static_cast<size_t>(end - start);
        double finish_percent = static_cast<double>(progress_counter)
                                / static_cast<double>(n_snps) * 100;

        pbar.status->message(
            fmt::format(
                "{:.1f}% ({}/{}) | {}",
                finish_percent,
                gelex::HumanReadable(n_tested),
                gelex::HumanReadable(n_snps),
                eta_calculator.get_eta(n_tested)));
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
