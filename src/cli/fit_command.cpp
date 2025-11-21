#include "gelex/cli/fit_command.h"

#include <thread>
#include "fit_command_detail.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/data_pipe.h"
#include "gelex/data/sample_manager.h"
#include "gelex/estimator/bayes/params.h"
#include "gelex/logger.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/effects.h"

void fit_command(argparse::ArgumentParser& cmd)
{
    cmd.add_description(
        "Fit genomic prediction models using Bayesian or GBLUP methods");

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
    cmd.add_argument("--qcovar")
        .default_value("")
        .help("Quantitative covariates (TSV: FID, IID, covar1, ...)");
    cmd.add_argument("--covar").default_value("").help(
        "Categorical covariates (TSV: FID, IID, factor1, ...)");
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
        .help("SNPs per chunk (controls memory usage)")
        .default_value(10000)
        .scan<'i', int>();
    cmd.add_argument("--iid-only")
        .help("Use only IID for sample matching (ignore FID)")
        .flag();

    // ================================================================
    // Model Configuration
    // ================================================================
    cmd.add_group("Model Configuration");
    cmd.add_argument("-m", "--method")
        .help(
            "Method: A/B/C/R/RR/GBLUP (+d for dominance, +pi to estimate "
            "mixture), e.g. RRd, Bdpi")
        .default_value("RR")
        .metavar("<METHOD>")
        .choices(
            "A",
            "Ad",
            "B",
            "Bpi",
            "Bd",
            "Bdpi",
            "C",
            "Cpi",
            "Cd",
            "Cdpi",
            "R",
            "Rd",
            "RR",
            "RRd",
            "GBLUP")
        .required();

    cmd.add_argument("--scale")
        .help("Additive variance scales for BayesR (5 values)")
        .nargs(5)
        .default_value(std::vector<double>{0, 0.001, 0.01, 0.1, 1})
        .scan<'g', double>();
    cmd.add_argument("--pi")
        .help("Additive mixture proportions for BayesB/C/R")
        .nargs(2)
        .default_value(std::vector<double>{0.95, 0.05})
        .scan<'g', double>();
    cmd.add_argument("--dscale")
        .help("Dominance variance scales for BayesR (5 values)")
        .nargs(5)
        .default_value(std::vector<double>{0, 0.001, 0.01, 0.1, 1})
        .scan<'g', double>();
    cmd.add_argument("--dpi")
        .help("Dominance mixture proportions for BayesB/C/R")
        .nargs(2)
        .default_value(std::vector<double>{0.95, 0.05})
        .scan<'g', double>();

    // ================================================================
    // MCMC Parameters
    // ================================================================
    cmd.add_group("MCMC Configuration");
    cmd.add_argument("--iters")
        .help("Total MCMC iterations")
        .default_value(5000)
        .scan<'i', int>();
    cmd.add_argument("--burnin")
        .help("Burn-in iterations to discard")
        .default_value(4000)
        .scan<'i', int>();
    cmd.add_argument("--thin")
        .help("Thinning interval for samples")
        .default_value(1)
        .scan<'i', int>();
    cmd.add_argument("--chains")
        .help("Number of MCMC chains")
        .default_value(1)
        .scan<'i', int>();

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
    cmd.add_argument("--mmap")
        .help(
            "Use memory-mapped I/O for genotype matrix(much lower RAM, may be "
            "slower)")
        .flag();
}

int fit_execute(argparse::ArgumentParser& fit)
{
    // ================================================================
    // ====================== Preparations ============================
    // ================================================================
    std::string out_prefix = fit.get("--out");
    gelex::BayesAlphabet type = gelex::get_bayesalphabet(fit.get("-m"))
                                    .value_or(gelex::BayesAlphabet::RR);
    bool dom = app::has_dominance(type);

    gelex::logging::initialize(out_prefix);
    auto logger = gelex::logging::get();

    app::setup_parallelization(fit.get<int>("--threads"), logger);

    auto bed_result = gelex::BedPipe::format_bed_path(fit.get("bfile"));
    VALIDATE_RESULT_OR_RETURN(bed_result, logger);
    auto bed_path = bed_result.value();

    gelex::DataPipe::Config config{
        .phenotype_path = fit.get("pheno"),
        .phenotype_column = fit.get<int>("--pheno-col"),
        .bed_path = bed_path,
        .use_dominance_effect = dom,
        .use_mmap = fit.get<bool>("--mmap"),
        .chunk_size = fit.get<int>("--chunk-size"),
        .qcovar_path = fit.get("--qcovar"),
        .covar_path = fit.get("--covar"),
        .iid_only = fit.get<bool>("--iid-only"),
        .output_prefix = fit.get("--out")};

    auto data_pipe = gelex::DataPipe::create(config);
    VALIDATE_RESULT_OR_RETURN(data_pipe, logger);

    auto model_result = gelex::BayesModel::create(*data_pipe);
    VALIDATE_RESULT_OR_RETURN(model_result, logger);
    auto model = std::move(model_result.value());

    if (app::configure_model_priors(model, type, fit, logger) != 0)
    {
        return 1;
    }

    const gelex::MCMCParams mcmc_params(
        fit.get<int>("--iters"),
        fit.get<int>("--burnin"),
        fit.get<int>("--thin"),
        fit.get<int>("--chains"));

    if (app::run_mcmc_analysis(
            model,
            type,
            mcmc_params,
            bed_path.replace_extension(".bim"),
            out_prefix,
            logger)
        != 0)
    {
        return 1;
    }

    return 0;
}
