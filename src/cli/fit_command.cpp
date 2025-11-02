#include "gelex/cli/fit_command.h"

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
        "Fit a genomic prediction model using Bayesian or GBLUP methods");
    cmd.add_argument("--pheno")
        .help("phenotype file with columns [FID, IID, A, B...], sep by tab")
        .required();
    cmd.add_argument("--qcovar")
        .default_value("")
        .help(
            "quantitative covariate file with columns [FID, IID, C1, C2...], "
            "sep by tab");
    cmd.add_argument("--covar").default_value("").help(
        "categorical covariate file with columns [FID, IID, C1, C2...], sep by "
        "tab");
    cmd.add_argument("--bfile")
        .help("prefix for PLINK1 binary files")
        .required();
    cmd.add_argument("--pheno-col")
        .help(
            "specify which phenotype column to use, default is the 3rd column")
        .default_value(3)
        .scan<'i', int>();
    cmd.add_argument("--chunk-size")
        .help("chunk size for processing snps, default is 10000")
        .default_value(10000)
        .scan<'i', int>();
    cmd.add_argument("--dr-mean")
        .help("prior mean for dominance to additive ratio, default is 0.0")
        .default_value(0.0)
        .scan<'g', double>();
    cmd.add_argument("--dr-var")
        .help("prior variance for dominance to additive ratio, default is 1.0")
        .default_value(1.0)
        .scan<'g', double>();
    cmd.add_argument("--scale")
        .help(
            "variance scales for additive mixture components (e.g., for "
            "BayesR: --scale "
            "0.0001 0.001 0.01 0.1 1.0)")
        .nargs(5)
        .scan<'g', double>();

    cmd.add_argument("--pi")
        .help("mixture proportions for additive effects (e.g., --pi 0.95 0.05)")
        .nargs(2)
        .scan<'g', double>();

    cmd.add_argument("--dscale")
        .help("variance scales for dominance mixture components")
        .nargs(5)
        .scan<'g', double>();

    cmd.add_argument("--dpi")
        .help(
            "mixture proportions for dominance effects (e.g., --dpi 0.95 0.05)")
        .nargs(2)
        .scan<'g', double>();
    cmd.add_argument("-m", "--method")
        .help(
            "genomic prediction method: A, Ad, B, Bpi, Bd, Bdpi, C, Cpi, Cd, "
            "Cdpi, R, Rd, RR, RRd, GBLUP")
        .default_value("RR")
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
    cmd.add_argument("-o", "--out").help("output prefix").required();
    cmd.add_argument("--iid_only")
        .help("use IID for sample ID, default is false, which will use FID_IID")
        .flag();
    cmd.add_argument("--iters")
        .help("number of total MCMC iterations, default is 3000")
        .default_value(3000)
        .scan<'i', int>();
    cmd.add_argument("--burnin")
        .help("number of burn-in iterations, default is 1000")
        .default_value(1000)
        .scan<'i', int>();
    cmd.add_argument("--thin")
        .help("thinning interval, default is 1")
        .default_value(1)
        .scan<'i', int>();
    cmd.add_argument("--chains")
        .help("number of MCMC chains, default is 1")
        .default_value(1)
        .scan<'i', int>();
    cmd.add_argument("--threads")
        .help("number of threads for parallelization, default is 16")
        .default_value(16)
        .scan<'i', int>();
    cmd.add_argument("--mmap")
        .help("use memory-mapped genotype files instead of in-memory loading")
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
    const bool use_mmap = fit.get<bool>("--mmap");
    const int chunk_size = fit.get<int>("--chunk-size");
    bool dom = app::has_dominance(type);

    gelex::logging::initialize(out_prefix);
    auto logger = gelex::logging::get();

    app::setup_parallelization(fit.get<int>("--threads"), logger);

    auto bed_result = gelex::valid_bed(fit.get("--bfile"));
    VALIDATE_RESULT_OR_RETURN(bed_result, logger);
    auto bed_path = bed_result.value();

    auto sample_manager_result = gelex::SampleManager::create(
        bed_path.replace_extension(".fam"), fit.get<bool>("--iid_only"));
    VALIDATE_RESULT_OR_RETURN(sample_manager_result, logger);

    auto sample_manager = std::make_shared<gelex::SampleManager>(
        std::move(*sample_manager_result));

    gelex::DataPipe::Config config{
        .phenotype_path = fit.get("--pheno"),
        .phenotype_column = fit.get<int>("--pheno-col"),
        .qcovar_path = fit.get("--qcovar"),
        .covar_path = fit.get("--covar"),
        .iid_only = fit.get<bool>("--iid_only"),
        .output_prefix = fit.get("--out")};

    auto data_pipe = gelex::DataPipe::create(config, sample_manager);
    VALIDATE_RESULT_OR_RETURN(data_pipe, logger);

    auto model_result = gelex::BayesModel::create(*data_pipe);
    VALIDATE_RESULT_OR_RETURN(model_result, logger);
    auto model = std::move(model_result.value());

    if (app::process_genotype_effect(
            model,
            sample_manager,
            bed_path.replace_extension(".bed"),
            out_prefix,
            chunk_size,
            use_mmap,
            /*is_dominance=*/false,
            logger)
        != 0)
    {
        logger->error("Failed to process additive genotype effect");
        return 1;
    }

    if (dom)
    {
        if (app::process_genotype_effect(
                model,
                sample_manager,
                bed_path.replace_extension(".bed"),
                out_prefix,
                chunk_size,
                use_mmap,
                /*is_dominance=*/true,
                logger)
            != 0)
        {
            logger->error("Failed to process dominance genotype effect");
            return 1;
        }
    }

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
