#include "gelex/cli/fit_command.h"

#include <omp.h>
#include <memory>

#include "gelex/data/data_pipe.h"
#include "gelex/data/genotype_loader.h"
#include "gelex/data/genotype_pipe.h"
#include "gelex/data/sample_manager.h"
#include "gelex/error.h"
#include "gelex/estimator/bayes/mcmc.h"
#include "gelex/estimator/bayes/params.h"
#include "gelex/estimator/bayes/result_writer.h"
#include "gelex/logger.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior_strategies.h"
#include "gelex/model/bayes/trait_model.h"
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
    bool dom = fit.get<bool>("--dom");
    gelex::logging::initialize(out_prefix);
    auto logger = gelex::logging::get();

    // ================================================================
    // ====================== Set Parallelization =====================
    // ================================================================
    int num_threads = fit.get<int>("--threads");
    if (num_threads > 0)
    {
        omp_set_num_threads(num_threads);
        Eigen::setNbThreads(num_threads);
        logger->info("Using {} threads for parallel computation", num_threads);
    }

    // ================================================================
    // ======================== Data Clean ============================
    // ================================================================
    auto bed = gelex::valid_bed(fit.get("--bfile"));
    if (!bed)
    {
        logger->error(bed.error().message);
        return 1;
    }

    // Sample Intersection
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

    // ================================================================
    // ======================== Model Container =======================
    // ================================================================

    auto model = gelex::BayesModel::create(*data_pipe);
    if (!model)
    {
        logger->error(model.error().message);
        return 1;
    }
    auto process_genotype_effect = [&](bool dominance) -> int
    {
        if (fit.get<bool>("--mmap"))
        {
            auto genotype_pipe = gelex::GenotypePipe::create(
                bed->replace_extension(".bed"),
                sample_manager,
                fit.get("--out"),
                dominance);
            if (!genotype_pipe
                && genotype_pipe.error().code
                       != gelex::ErrorCode::OutputFileExists)
            {
                logger->error(genotype_pipe.error().message);
                return 1;
            }
            if (!genotype_pipe
                && genotype_pipe.error().code
                       == gelex::ErrorCode::OutputFileExists)
            {
                if (dominance)
                {
                    auto result
                        = gelex::GenotypeMap::create(out_prefix + ".dom.bmat");
                    if (!result)
                    {
                        logger->error(result.error().message);
                        return 1;
                    }
                    (*model).add_dominance(std::move(result.value()));
                }
                if (!dominance)
                {
                    auto result
                        = gelex::GenotypeMap::create(out_prefix + ".add.bmat");
                    if (!result)
                    {
                        logger->error(result.error().message);
                        return 1;
                    }
                    (*model).add_additive(std::move(result.value()));
                }
            }

            if (genotype_pipe)
            {
                auto process_result
                    = genotype_pipe->process(fit.get<int>("--chunk-size"));
                if (!process_result)
                {
                    logger->error(process_result.error().message);
                    return 1;
                }
                if (dominance)
                {
                    (*model).add_dominance(std::move(process_result.value()));
                }
                else
                {
                    (*model).add_additive(std::move(process_result.value()));
                }
            }
        }
        else
        {
            auto genotype_loader = gelex::GenotypeLoader::create(
                bed->replace_extension(".bed"), sample_manager, dominance);
            if (!genotype_loader)
            {
                logger->error(genotype_loader.error().message);
                return 1;
            }
            auto process_result
                = genotype_loader->process(fit.get<int>("--chunk-size"));

            if (dominance)
            {
                (*model).add_dominance(std::move(process_result.value()));
            }
            else
            {
                (*model).add_additive(std::move(process_result.value()));
            }
        }
        return 0;
    };

    if (process_genotype_effect(false) != 0)
    {
        logger->error("Failed to process additive genotype effect");
        return 1;
    }

    if (dom)
    {
        if (process_genotype_effect(true) != 0)
        {
            logger->error("Failed to process dominance genotype effect");
            return 1;
        }
    }

    // ================================================================
    // ======================== Prior Management ======================
    // ================================================================
    auto prior_strategy = gelex::create_prior_strategy(type);
    if (!prior_strategy)
    {
        logger->error(
            "Failed to create prior strategy for model type: {}",
            fit.get("-m"));
        return 1;
    }

    // Create PriorConfig from command-line arguments
    gelex::PriorConfig prior_config;

    // Set default phenotype variance (this should be computed from actual data)
    // For now, we'll use a placeholder - this should be computed from phenotype
    // data
    prior_config.phenotype_variance = model->phenotype_variance();

    // Handle additive effect parameters
    if (auto pi_args = fit.present<std::vector<double>>("--pi"))
    {
        prior_config.additive.mixture_proportions
            = Eigen::Map<const Eigen::VectorXd>(
                pi_args->data(), pi_args->size());
    }
    if (auto scale_args = fit.present<std::vector<double>>("--scale"))
    {
        prior_config.additive.mixture_scales
            = Eigen::Map<const Eigen::VectorXd>(
                scale_args->data(), scale_args->size());
    }

    // Handle dominance effect parameters
    if (auto dpi_args = fit.present<std::vector<double>>("--dpi"))
    {
        prior_config.dominant.mixture_proportions
            = Eigen::Map<const Eigen::VectorXd>(
                dpi_args->data(), dpi_args->size());
    }
    if (auto dscale_args = fit.present<std::vector<double>>("--dscale"))
    {
        prior_config.dominant.mixture_scales
            = Eigen::Map<const Eigen::VectorXd>(
                dscale_args->data(), dscale_args->size());
    }

    // Set default mixture proportions and scales based on model type
    // Only set defaults if user didn't provide explicit values
    if (!fit.present<std::vector<double>>("--pi"))
    {
        switch (type)
        {
            using bt = gelex::BayesAlphabet;
            case bt::B:
            case bt::Bpi:
            case bt::Bd:
            case bt::Bdpi:
            case bt::C:
            case bt::Cpi:
            case bt::Cd:
            case bt::Cdpi:
                // Default for B/C models: {0.95, 0.05}
                prior_config.additive.mixture_proportions
                    = Eigen::Vector2d(0.95, 0.05);
                break;
            case bt::R:
            case bt::Rd:
                // Default for R models: {0.95, 0.02, 0.01, 0.01, 0.01}
                prior_config.additive.mixture_proportions = Eigen::VectorXd(5);
                prior_config.additive.mixture_proportions << 0.95, 0.02, 0.01,
                    0.01, 0.01;
                break;
            default:
                // For non-mixture models (A, RR, etc.), don't set pi
                break;
        }
    }

    // Set default scales for R models if not provided
    if (!fit.present<std::vector<double>>("--scale"))
    {
        switch (type)
        {
            using bt = gelex::BayesAlphabet;
            case bt::R:
            case bt::Rd:
                // Default scales for R models: {0, 0.001, 0.01, 0.1, 1}
                prior_config.additive.mixture_scales = Eigen::VectorXd(5);
                prior_config.additive.mixture_scales << 0.0, 0.001, 0.01, 0.1,
                    1.0;
                break;
            default:
                // For non-R models, don't set scales
                break;
        }
    }

    // Set default dominance mixture proportions if not provided
    if (!fit.present<std::vector<double>>("--dpi"))
    {
        switch (type)
        {
            using bt = gelex::BayesAlphabet;
            case bt::Bd:
            case bt::Bdpi:
            case bt::Cd:
            case bt::Cdpi:
                // Default for dominance B/C models: {0.95, 0.05}
                prior_config.dominant.mixture_proportions
                    = Eigen::Vector2d(0.95, 0.05);
                break;
            case bt::Rd:
                // Default for dominance R models: {0.95, 0.02, 0.01, 0.01,
                // 0.01}
                prior_config.dominant.mixture_proportions = Eigen::VectorXd(5);
                prior_config.dominant.mixture_proportions << 0.95, 0.02, 0.01,
                    0.01, 0.01;
                break;
            default:
                // For non-dominance models, don't set dpi
                break;
        }
    }

    // Set default dominance scales for Rd models if not provided
    if (!fit.present<std::vector<double>>("--dscale"))
    {
        switch (type)
        {
            using bt = gelex::BayesAlphabet;
            case bt::Rd:
                // Default scales for Rd models: {0, 0.001, 0.01, 0.1, 1}
                prior_config.dominant.mixture_scales = Eigen::VectorXd(5);
                prior_config.dominant.mixture_scales << 0.0, 0.001, 0.01, 0.1,
                    1.0;
                break;
            default:
                // For non-Rd models, don't set dscale
                break;
        }
    }

    // Apply prior configuration
    auto prior_result = (*prior_strategy)(*model, prior_config);
    if (!prior_result)
    {
        logger->error(
            "Failed to apply prior configuration: {}",
            prior_result.error().message);
        return 1;
    }

    // ================================================================
    // ========================== MCMC ================================
    // ================================================================
    gelex::MCMCParams mcmc_params(
        fit.get<int>("--iters"),
        fit.get<int>("--burnin"),
        fit.get<int>("--thin"),
        fit.get<int>("--chains"));

    auto bim_path = bed->replace_extension(".bim");
    auto run_and_write = [&](auto trait_model)
    {
        gelex::MCMC mcmc(mcmc_params, trait_model);
        gelex::MCMCResult result = mcmc.run(model.value());
        gelex::MCMCResultWriter writer(result, bim_path);
        writer.save(out_prefix);
    };

    // Select appropriate TraitModel based on method
    switch (type)
    {
        using bt = gelex::BayesAlphabet;
        case (bt::A):
            run_and_write(gelex::BayesA{});
            break;
        case (bt::Ad):
            run_and_write(gelex::BayesAd{});
            break;
        case (bt::B):
            run_and_write(gelex::BayesB{});
            break;
        case (bt::Bpi):
            run_and_write(gelex::BayesBpi{});
            break;
        case (bt::Bd):
            run_and_write(gelex::BayesBD{});
            break;
        case (bt::Bdpi):
            run_and_write(gelex::BayesBdpi{});
            break;
        case (bt::C):
            run_and_write(gelex::BayesC{});
            break;
        case (bt::Cpi):
            run_and_write(gelex::BayesCpi{});
            break;
        case (bt::Cd):
            run_and_write(gelex::BayesCd{});
            break;
        case (bt::Cdpi):
            run_and_write(gelex::BayesCdpi{});
            break;
        case (bt::R):
            run_and_write(gelex::BayesR{});
            break;
        case (bt::Rd):
            run_and_write(gelex::BayesRd{});
            break;
        case (bt::RR):
            run_and_write(gelex::BayesRR{});
            break;
        case (bt::RRd):
            run_and_write(gelex::BayesRRD{});
            break;
        default:
            logger->error("Unsupported method: {}", fit.get("-m"));
            break;
    }

    return 0;
}
