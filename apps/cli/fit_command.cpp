#include "fit_command.h"

#include <barkeep.h>
#include <fmt/format.h>
#include <omp.h>
#include <spdlog/logger.h>

#include <Eigen/Core>
#include <filesystem>
#include <memory>
#include <thread>
#include <vector>

#include "cli_helper.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/data_pipe.h"
#include "gelex/estimator/bayes/mcmc.h"
#include "gelex/estimator/bayes/params.h"
#include "gelex/estimator/bayes/result_writer.h"
#include "gelex/logger.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior_strategies.h"
#include "gelex/model/bayes/trait_model.h"
#include "gelex/model/effects.h"
#include "utils/formatter.h"

namespace bk = barkeep;

namespace
{

using gelex::BayesAlphabet;

auto has_dominance(BayesAlphabet type) -> bool
{
    switch (type)
    {
        case BayesAlphabet::Bd:
        case BayesAlphabet::Bdpi:
        case BayesAlphabet::Cd:
        case BayesAlphabet::Cdpi:
        case BayesAlphabet::Rd:
        case BayesAlphabet::Ad:
        case BayesAlphabet::RRd:
            return true;
        default:
            return false;
    }
}

auto get_default_pi(BayesAlphabet type) -> Eigen::VectorXd
{
    switch (type)
    {
        case BayesAlphabet::B:
        case BayesAlphabet::Bpi:
        case BayesAlphabet::Bd:
        case BayesAlphabet::Bdpi:
        case BayesAlphabet::C:
        case BayesAlphabet::Cpi:
        case BayesAlphabet::Cd:
        case BayesAlphabet::Cdpi:
            return Eigen::Vector2d(0.95, 0.05);
        case BayesAlphabet::R:
        case BayesAlphabet::Rd:
        {
            Eigen::VectorXd v(5);
            v << 0.95, 0.02, 0.01, 0.01, 0.01;
            return v;
        }
        case BayesAlphabet::A:
        case BayesAlphabet::RR:
        case BayesAlphabet::Ad:
        case BayesAlphabet::RRd:
            return Eigen::Vector2d(0.0, 1.0);
        default:
            return Eigen::VectorXd{};
    }
}

auto get_default_scale(BayesAlphabet type) -> Eigen::VectorXd
{
    switch (type)
    {
        case BayesAlphabet::R:
        case BayesAlphabet::Rd:
        {
            Eigen::VectorXd v(5);
            v << 0.0, 0.001, 0.01, 0.1, 1.0;
            return v;
        }
        default:
            return Eigen::VectorXd{};
    }
}

auto configure_model_priors(
    gelex::BayesModel& model,
    BayesAlphabet type,
    argparse::ArgumentParser& fit,
    const std::shared_ptr<spdlog::logger>& logger) -> int
{
    auto prior_strategy = gelex::create_prior_strategy(type);
    if (!prior_strategy)
    {
        logger->error(
            "Failed to create prior strategy for model type: {}",
            fit.get("-m"));
        return 1;
    }

    gelex::PriorConfig prior_config;
    prior_config.phenotype_variance = model.phenotype_variance();

    // Configure Additive Pi
    if (fit.is_used("--pi"))
    {
        auto pi = fit.get<std::vector<double>>("--pi");
        prior_config.additive.mixture_proportions
            = Eigen::Map<const Eigen::VectorXd>(
                pi.data(), static_cast<Eigen::Index>(pi.size()));
    }
    else
    {
        prior_config.additive.mixture_proportions = get_default_pi(type);
    }

    // Configure Dominance Pi
    if (fit.is_used("--dpi"))
    {
        auto dpi = fit.get<std::vector<double>>("--dpi");
        prior_config.dominant.mixture_proportions
            = Eigen::Map<const Eigen::VectorXd>(
                dpi.data(), static_cast<Eigen::Index>(dpi.size()));
    }
    else
    {
        prior_config.dominant.mixture_proportions = get_default_pi(type);
    }

    // Configure Additive Scale
    if (fit.is_used("--scale"))
    {
        auto scale = fit.get<std::vector<double>>("--scale");
        prior_config.additive.mixture_scales
            = Eigen::Map<const Eigen::VectorXd>(
                scale.data(), static_cast<Eigen::Index>(scale.size()));
    }
    else
    {
        prior_config.additive.mixture_scales = get_default_scale(type);
    }

    // Configure Dominance Scale
    if (fit.is_used("--dscale"))
    {
        auto dscale = fit.get<std::vector<double>>("--dscale");
        prior_config.dominant.mixture_scales
            = Eigen::Map<const Eigen::VectorXd>(
                dscale.data(), static_cast<Eigen::Index>(dscale.size()));
    }
    else
    {
        prior_config.dominant.mixture_scales = get_default_scale(type);
    }

    (*prior_strategy)(model, prior_config);
    return 0;
}

auto run_mcmc_analysis(
    gelex::BayesModel& model,
    BayesAlphabet type,
    const gelex::MCMCParams& mcmc_params,
    const std::filesystem::path& bim_path,
    std::string_view out_prefix,
    const std::shared_ptr<spdlog::logger>& logger) -> int
{
    auto run_and_write = [&](auto trait_model)
    {
        gelex::MCMC mcmc(mcmc_params, trait_model);
        gelex::MCMCResult result = mcmc.run(model);
        gelex::MCMCResultWriter writer(result, bim_path);
        writer.save(out_prefix);
    };

    switch (type)
    {
        case (BayesAlphabet::A):
            run_and_write(gelex::BayesA{});
            break;
        case (BayesAlphabet::Ad):
            run_and_write(gelex::BayesAd{});
            break;
        case (BayesAlphabet::B):
            run_and_write(gelex::BayesB{});
            break;
        case (BayesAlphabet::Bpi):
            run_and_write(gelex::BayesBpi{});
            break;
        case (BayesAlphabet::Bd):
            run_and_write(gelex::BayesBd{});
            break;
        case (BayesAlphabet::Bdpi):
            run_and_write(gelex::BayesBdpi{});
            break;
        case (BayesAlphabet::C):
            run_and_write(gelex::BayesC{});
            break;
        case (BayesAlphabet::Cpi):
            run_and_write(gelex::BayesCpi{});
            break;
        case (BayesAlphabet::Cd):
            run_and_write(gelex::BayesCd{});
            break;
        case (BayesAlphabet::Cdpi):
            run_and_write(gelex::BayesCdpi{});
            break;
        case (BayesAlphabet::R):
            run_and_write(gelex::BayesR{});
            break;
        case (BayesAlphabet::Rd):
            run_and_write(gelex::BayesRd{});
            break;
        case (BayesAlphabet::RR):
            run_and_write(gelex::BayesRR{});
            break;
        case (BayesAlphabet::RRd):
            run_and_write(gelex::BayesRRd{});
            break;
        default:
            logger->error("Unsupported method: {}", type);
            return 1;
    }
    return 0;
}

}  // namespace

int fit_execute(argparse::ArgumentParser& fit)
{
    // ================================================================
    // ====================== Preparations ============================
    // ================================================================
    std::string out_prefix = fit.get("--out");
    gelex::BayesAlphabet type = gelex::get_bayesalphabet(fit.get("-m"))
                                    .value_or(gelex::BayesAlphabet::RR);
    bool dom = has_dominance(type);

    auto logger = gelex::logging::get();

    gelex::cli::setup_parallelization(fit.get<int>("--threads"));

    gelex::cli::print_fit_header(
        fit.get("-m"),
        dom,
        fit.get<int>("--iters"),
        fit.get<int>("--burnin"),
        fit.get<int>("--threads"));

    auto bed_path = gelex::BedPipe::format_bed_path(fit.get("bfile"));

    gelex::DataPipe::Config config{
        .phenotype_path = fit.get("pheno"),
        .phenotype_column = fit.get<int>("--pheno-col"),
        .bed_path = bed_path,
        .use_dominance_effect = dom,
        .use_mmap = fit.get<bool>("--mmap"),
        .chunk_size = fit.get<int>("--chunk-size"),
        .qcovar_path = fit.get("--qcovar"),
        .dcovar_path = fit.get("--dcovar"),
        .iid_only = fit.get<bool>("--iid-only"),
        .output_prefix = fit.get("--out")};

    // ================================================================
    // Data Loading & Pipeline
    // ================================================================
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
            "No common samples found between phenotype, covariates, and "
            "genotype files.");
        return 1;
    }

    logger->info(gelex::task("Matrix Construction:"));
    logger->info(gelex::subtask("Additive:"));
    auto add_stats = data_pipe.load_additive_matrix();

    logger->info(gelex::subsubtask("{} SNPs processed", add_stats.num_snps));
    logger->info(
        gelex::subsubtask(
            "{} monomorphic SNPs excluded", add_stats.monomorphic_snps));

    if (config.use_dominance_effect)
    {
        logger->info(gelex::subtask("Dominance:"));
        auto dom_stats = data_pipe.load_dominance_matrix();

        logger->info(
            gelex::subsubtask("{} SNPs processed", dom_stats.num_snps));
        logger->info(
            gelex::subsubtask(
                "{} monomorphic SNPs excluded", dom_stats.monomorphic_snps));
    }

    data_pipe.finalize();

    gelex::BayesModel model(data_pipe);

    if (configure_model_priors(model, type, fit, logger) != 0)
    {
        return 1;
    }

    const gelex::MCMCParams mcmc_params(
        fit.get<int>("--iters"),
        fit.get<int>("--burnin"),
        fit.get<int>("--thin"),
        fit.get<int>("--chains"));

    if (run_mcmc_analysis(
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
    logger->info(gelex::success("Parameters saved to  : {}.param", out_prefix));
    logger->info(
        gelex::success("SNP Effects saved to : {}.snp.eff", out_prefix));
    logger->info(gelex::success("Run Log saved to     : {}.log", out_prefix));
    logger->info(
        fmt::format(
            fmt::emphasis::bold | fmt::fg(fmt::color::light_cyan),
            "───────────────────────────────────"
            "───────────────────────────────────"));

    return 0;
}
