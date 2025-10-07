#include "gelex/cli/fit_command.h"
#include "gelex/data/data_pipe.h"
#include "gelex/data/genotype_pipe.h"
#include "gelex/data/sample_manager.h"
#include "gelex/error.h"
#include "gelex/estimator/bayes/mcmc.h"
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
    cmd.add_argument("--dom")
        .help("enable estimation of dominance effects")
        .flag();
    cmd.add_argument("-m", "--method")
        .help("genomic prediction method: A, B(pi), C(pi), R, RR, GBLUP")
        .default_value("RR")
        .choices("A", "B", "Bpi", "C", "Cpi", "R", "RR", "GBLUP")
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
}

int fit_excute(argparse::ArgumentParser& fit)
{
    std::string out_prefix = fit.get("--out");
    gelex::BayesAlphabet type = gelex::get_bayesalphabet(fit.get("-m"))
                                    .value_or(gelex::BayesAlphabet::RR);

    gelex::logging::initialize(out_prefix);
    auto logger = gelex::logging::get();

    auto bed = gelex::valid_bed(fit.get("--bfile"));

    if (!bed)
    {
        logger->error(bed.error().message);
        return 1;
    }

    // Create shared SampleManager
    auto sample_manager_result = gelex::SampleManager::create(
        bed->replace_extension(".fam"), fit.get<bool>("--iid_only"));

    if (!sample_manager_result)
    {
        logger->error(sample_manager_result.error().message);
        return 1;
    }

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

    if (!data_pipe)
    {
        logger->error(data_pipe.error().message);
        return 1;
    }

    auto model = gelex::BayesModel::create_from_datapipe(*data_pipe, type);
    {
        auto genotype_pipe = gelex::GenotypePipe::create(
            bed->replace_extension(".bed"), sample_manager, fit.get("--out"));

        if (!genotype_pipe
            && genotype_pipe.error().code != gelex::ErrorCode::OutputFileExists)
        {
            logger->error(genotype_pipe.error().message);
            return 1;
        }

        if (genotype_pipe)
        {
            if (auto process_result
                = genotype_pipe->process<gelex::NonStandardizingProcessor>(
                    fit.get<int>("--chunk-size"));
                !process_result)

            {
                logger->error(process_result.error().message);
                return 1;
            }
        }

        auto gmap = gelex::GenotypeMap::create(out_prefix + ".add.bmat");
        if (!gmap)
        {
            logger->error(gmap.error().message);
            return 1;
        }
        model->add_additive_effect(std::move(gmap.value()));
    }

    if (fit.get<bool>("--dom"))
    {
        auto genotype_pipe = gelex::GenotypePipe::create(
            bed->replace_extension(".bed"),
            sample_manager,
            fit.get("--out"),
            true);
        if (!genotype_pipe
            && genotype_pipe.error().code != gelex::ErrorCode::OutputFileExists)
        {
            logger->error(genotype_pipe.error().message);
            return 1;
        }

        if (genotype_pipe)
        {
            if (auto process_result
                = genotype_pipe->process<gelex::NonStandardizingProcessor>(
                    fit.get<int>("--chunk-size"));
                !process_result)

            {
                logger->error(process_result.error().message);
                return 1;
            }
        }

        auto gmap = gelex::GenotypeMap::create(out_prefix + ".dom.bmat");
        if (!gmap)
        {
            logger->error(gmap.error().message);
            return 1;
        }
        model->add_dominance_effect(std::move(gmap.value()));
    }

    gelex::MCMCParams mcmc_params(
        fit.get<int>("--iters"),
        fit.get<int>("--burnin"),
        fit.get<int>("--thin"),
        fit.get<int>("--chains"));

    gelex::MCMC mcmc(mcmc_params);
    mcmc.run(model.value(), out_prefix);
    return 0;
}
