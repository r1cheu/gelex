#include "gelex/cli/gwas_command.h"

#include <fmt/format.h>
#include <omp.h>
#include <memory>
#include <thread>

#include "../src/utils/formatter.h"
#include "config.h"
#include "gelex/cli/utils.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/data_pipe.h"
#include "gelex/estimator/freq/estimator.h"
#include "gelex/gwas/association_test.h"
#include "gelex/gwas/gwas_writer.h"
#include "gelex/gwas/snp_encoder.h"
#include "gelex/logger.h"
#include "gelex/model/freq/model.h"
#include "gelex/optim/optimizer_state.h"
#include "gelex/optim/variance_calculator.h"

#include "../src/data/loader/bim_loader.h"

namespace
{

void setup_parallelization(int num_threads)
{
    if (num_threads > 0)
    {
        omp_set_num_threads(num_threads);
        Eigen::setNbThreads(num_threads);
    }
}

}  // namespace

auto gwas_command(argparse::ArgumentParser& cmd) -> void
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
    // Data Processing
    // ================================================================
    cmd.add_group("Processing Options");
    cmd.add_argument("--pheno-col")
        .help("Phenotype column index (0-based)")
        .default_value(2)
        .scan<'i', int>();
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
    cmd.add_argument("--test")
        .help("Test type for ad model: joint or separate")
        .default_value("joint")
        .metavar("<TEST>")
        .choices("joint", "separate");

    // ================================================================
    // Performance
    // ================================================================
    cmd.add_group("Performance");
    const int default_threads
        = std::max(1u, std::thread::hardware_concurrency() / 2);
    cmd.add_argument("--threads")
        .help("Number of CPU threads to use")
        .default_value(default_threads)
        .scan<'i', int>();
}

auto gwas_execute(argparse::ArgumentParser& cmd) -> int
{
    auto logger = gelex::logging::get();
    std::string out_prefix = cmd.get("--out");

    // Parse model and test type
    auto gwas_model = gelex::gwas::parse_gwas_model(cmd.get("--model"));
    auto test_type = gelex::gwas::parse_test_type(cmd.get("--test"));

    setup_parallelization(cmd.get<int>("--threads"));

    // Get GRM paths
    auto grm_paths = cmd.get<std::vector<std::string>>("--grm");

    logger->info(
        fmt::format(
            fmt::emphasis::bold | fmt::fg(fmt::color::light_cyan),
            "───────────────────────────────────"
            "───────────────────────────────────"));
    logger->info(
        fmt::format(
            fmt::emphasis::bold | fmt::fg(fmt::color::light_cyan),
            "                    GWAS Analysis - gelex v{}",
            PROJECT_VERSION));
    logger->info(
        fmt::format(
            fmt::emphasis::bold | fmt::fg(fmt::color::light_cyan),
            "───────────────────────────────────"
            "───────────────────────────────────"));

    logger->info(gelex::task("Model     : {}", cmd.get("--model")));
    if (gwas_model == gelex::gwas::GwasModel::AdditiveDominance)
    {
        logger->info(gelex::task("Test type : {}", cmd.get("--test")));
    }
    logger->info(gelex::task("Threads   : {}", cmd.get<int>("--threads")));

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
        .additive_grm_path = grm_paths[0],
        .dominance_grm_path = grm_paths.size() > 1 ? grm_paths[1] : ""};

    gelex::DataPipe data_pipe(config);

    logger->info("");
    logger->info(gelex::section("Loading Data..."));

    auto p_stats = data_pipe.load_phenotypes();
    logger->info(
        gelex::success(
            "Phenotypes : {} samples ('{}')",
            p_stats.samples_loaded,
            p_stats.trait_name));

    logger->info(
        gelex::success(
            "Genotypes  : {} samples", data_pipe.num_genotype_samples()));

    auto c_stats = data_pipe.load_covariates();
    if (c_stats.qcovar_loaded > 0 || c_stats.dcovar_loaded > 0)
    {
        logger->info(gelex::task("Covariates : "));
    }
    if (c_stats.qcovar_loaded > 0)
    {
        logger->info(
            gelex::subtask(
                "Quantitative : {} loaded ",
                gelex::format_names(c_stats.q_names)));
    }
    if (c_stats.dcovar_loaded > 0)
    {
        logger->info(
            gelex::subtask(
                "Discrete     : {} loaded ",
                gelex::format_names(c_stats.d_names)));
    }

    // Load GRM(s)
    logger->info(gelex::task("GRM        : "));
    auto grm_stats = data_pipe.load_additive_grm();
    logger->info(
        gelex::subtask("Additive   : {} samples", grm_stats.samples_in_file));

    if (grm_paths.size() > 1)
    {
        auto dom_grm_stats = data_pipe.load_dominance_grm();
        logger->info(
            gelex::subtask(
                "Dominance  : {} samples", dom_grm_stats.samples_in_file));
    }

    // Intersect samples
    logger->info("");
    logger->info(gelex::section("Pre-processing..."));
    auto i_stats = data_pipe.intersect_samples();
    logger->info(gelex::task("Sample Intersection:"));
    logger->info(
        gelex::subtask("Common samples : {} ", i_stats.common_samples));
    logger->info(
        gelex::subtask("Excluded       : {} ", i_stats.excluded_samples));

    if (i_stats.common_samples == 0)
    {
        logger->error(
            "No common samples found between phenotype, covariates, GRM, and "
            "genotype files.");
        return 1;
    }

    data_pipe.finalize();

    // ================================================================
    // Build null model and fit with REML
    // ================================================================
    logger->info("");
    logger->info(gelex::section("Fitting null model (REML)..."));

    gelex::FreqModel model(data_pipe);
    gelex::FreqState state(model);

    // Initialize variance components
    double pheno_var = model.phenotype().array().square().mean();
    state.residual().variance = pheno_var * 0.5;
    for (auto& g : state.genetic())
    {
        g.variance
            = pheno_var * 0.5 / static_cast<double>(state.genetic().size());
    }

    gelex::Estimator estimator(100, 1e-6);
    estimator.fit(model, state, true, true);

    if (!estimator.is_converged())
    {
        logger->warn("REML did not converge, results may be unreliable");
    }

    // Get V^{-1} from optimizer state
    // Re-compute it since EstimatorNew doesn't expose internal state
    gelex::OptimizerState opt_state(model);
    gelex::variance_calculator::compute_v(model, state, opt_state.v);
    gelex::variance_calculator::v_inv_logdet(opt_state.v);
    // Now opt_state.v contains V^{-1}

    // Compute residuals: r = y - Xβ̂
    const Eigen::VectorXd& y = model.phenotype();
    const Eigen::MatrixXd& X = model.fixed().X;
    const Eigen::VectorXd& beta = state.fixed().coeff;
    Eigen::VectorXd residual = y - X * beta;

    // ================================================================
    // Load BIM file for SNP info
    // ================================================================
    auto bim_path = bed_path;
    bim_path.replace_extension(".bim");
    gelex::detail::BimLoader bim_loader(bim_path);

    // Create BedPipe for chunk loading using the same sample manager
    auto fam_path = bed_path;
    fam_path.replace_extension(".fam");
    auto sample_manager
        = std::make_shared<gelex::SampleManager>(fam_path, config.iid_only);

    // Note: BedPipe will use all samples from fam file
    // In a full implementation, we should use the same sample ordering
    // as the null model. For now, we assume samples are properly aligned.
    gelex::BedPipe bed_pipe(bed_path, sample_manager);

    const Eigen::Index n_snps = bed_pipe.num_snps();
    const Eigen::Index n_samples
        = static_cast<Eigen::Index>(i_stats.common_samples);
    const int chunk_size = cmd.get<int>("--chunk-size");

    logger->info("");
    logger->info(gelex::section("Running association tests..."));
    logger->info(gelex::task("SNPs to test : {}", n_snps));
    logger->info(gelex::task("Chunk size   : {}", chunk_size));

    // ================================================================
    // Association testing
    // ================================================================
    gelex::gwas::GwasWriter writer(out_prefix, gwas_model, test_type);
    writer.write_header();

    size_t n_tested = 0;
    size_t n_chunks
        = static_cast<size_t>((n_snps + chunk_size - 1) / chunk_size);

    for (size_t chunk_idx = 0; chunk_idx < n_chunks; ++chunk_idx)
    {
        Eigen::Index start = static_cast<Eigen::Index>(chunk_idx)
                             * static_cast<Eigen::Index>(chunk_size);
        Eigen::Index end
            = std::min(start + static_cast<Eigen::Index>(chunk_size), n_snps);

        Eigen::MatrixXd geno_chunk = bed_pipe.load_chunk(start, end);

        // Store results for this chunk (to avoid critical section overhead)
        std::vector<
            std::pair<gelex::gwas::SNPInfo, gelex::gwas::AssociationResult>>
            chunk_results(static_cast<size_t>(end - start));

// Process each SNP in the chunk
#pragma omp parallel for schedule(dynamic)
        for (Eigen::Index snp_idx = 0; snp_idx < geno_chunk.cols(); ++snp_idx)
        {
            Eigen::Index global_idx = start + snp_idx;

            // Get raw genotype
            Eigen::VectorXd raw_geno = geno_chunk.col(snp_idx);

            // Encode SNP
            auto encoded = gelex::gwas::encode_snp(raw_geno, gwas_model);

            // Perform Wald test
            auto result = gelex::gwas::wald_test(
                encoded, residual, opt_state.v, gwas_model, test_type);

            // Get SNP info
            const auto& snp_meta
                = bim_loader.info()[static_cast<size_t>(global_idx)];

            gelex::gwas::SNPInfo snp_info{
                .chrom = snp_meta.chrom,
                .rsid = snp_meta.id,
                .bp = snp_meta.pos,
                .a1 = std::string(1, snp_meta.A1),
                .a2 = std::string(1, snp_meta.A2),
                .freq = encoded.maf,
                .n = static_cast<int>(n_samples)};

            chunk_results[static_cast<size_t>(snp_idx)]
                = std::make_pair(std::move(snp_info), result);
        }

        // Write results sequentially to maintain order
        for (const auto& [snp_info, result] : chunk_results)
        {
            writer.write_result(snp_info, result);
        }

        n_tested += static_cast<size_t>(end - start);

        if ((chunk_idx + 1) % 10 == 0 || chunk_idx == n_chunks - 1)
        {
            logger->info(
                gelex::subtask(
                    "Progress: {}/{} SNPs ({:.1f}%)",
                    n_tested,
                    n_snps,
                    100.0 * static_cast<double>(n_tested)
                        / static_cast<double>(n_snps)));
        }
    }

    writer.finalize();

    logger->info("");
    logger->info(gelex::success("GWAS complete!"));
    logger->info(gelex::success("Results saved to : {}.gwas.tsv", out_prefix));
    logger->info(
        fmt::format(
            fmt::emphasis::bold | fmt::fg(fmt::color::light_cyan),
            "───────────────────────────────────"
            "───────────────────────────────────"));

    return 0;
}
