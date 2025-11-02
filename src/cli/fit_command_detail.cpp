#include "fit_command_detail.h"

#include <omp.h>
#include <spdlog/logger.h>
#include <Eigen/Core>
#include "gelex/argparse.h"

#include "gelex/data/genotype_loader.h"
#include "gelex/data/genotype_matrix.h"
#include "gelex/data/genotype_mmap.h"
#include "gelex/data/genotype_pipe.h"
#include "gelex/estimator/bayes/mcmc.h"
#include "gelex/estimator/bayes/params.h"
#include "gelex/estimator/bayes/result_writer.h"
#include "gelex/model/bayes/prior_strategies.h"
#include "gelex/model/bayes/trait_model.h"

namespace app
{
auto has_dominance(gelex::BayesAlphabet type) -> bool
{
    switch (type)
    {
        using bt = gelex::BayesAlphabet;
        case bt::B:
        case bt::Bpi:
        case bt::C:
        case bt::Cpi:
        case bt::R:
        case bt::A:
        case bt::RR:
            return false;
        case bt::Bd:
        case bt::Bdpi:
        case bt::Cd:
        case bt::Cdpi:
        case bt::Rd:
        case bt::Ad:
        case bt::RRd:
            return true;
    }
}

auto set_default_pi_prior(gelex::BayesAlphabet type) -> Eigen::VectorXd
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
            return Eigen::VectorXd{{0.95, 0.05}};
        case bt::R:
        case bt::Rd:
            return Eigen::VectorXd{{0.95, 0.02, 0.01, 0.01, 0.01}};
        case gelex::BayesAlphabet::A:
        case gelex::BayesAlphabet::RR:
        case gelex::BayesAlphabet::Ad:
        case gelex::BayesAlphabet::RRd:
            return Eigen::VectorXd{{0.0, 1.0}};
    }
}

auto set_pi_prior(
    gelex::BayesAlphabet type,
    argparse::ArgumentParser& fit,
    gelex::PriorConfig& prior_config) -> void
{
    if (auto pi = fit.present<std::vector<double>>("--pi"))
    {
        prior_config.additive.mixture_proportions
            = Eigen::Map<const Eigen::VectorXd>(
                pi->data(), static_cast<Eigen::Index>(pi->size()));
    }
    else
    {
        prior_config.additive.mixture_proportions = set_default_pi_prior(type);
    }

    if (auto dpi = fit.present<std::vector<double>>("--dpi"))
    {
        prior_config.dominant.mixture_proportions
            = Eigen::Map<const Eigen::VectorXd>(
                dpi->data(), static_cast<Eigen::Index>(dpi->size()));
    }
    else
    {
        prior_config.dominant.mixture_proportions = set_default_pi_prior(type);
    }
}

auto set_default_scale_prior(gelex::BayesAlphabet type) -> Eigen::VectorXd
{
    switch (type)
    {
        using bt = gelex::BayesAlphabet;
        case bt::R:
        case bt::Rd:
            return Eigen::VectorXd{{0.0, 0.001, 0.01, 0.1, 1.0}};
        default:
            return Eigen::VectorXd{};
    }
}

auto set_scale_prior(
    gelex::BayesAlphabet type,
    argparse::ArgumentParser& fit,
    gelex::PriorConfig& prior_config) -> void
{
    if (auto scale = fit.present<std::vector<double>>("--scale"))
    {
        prior_config.additive.mixture_scales
            = Eigen::Map<const Eigen::VectorXd>(
                scale->data(), static_cast<Eigen::Index>(scale->size()));
    }
    else
    {
        prior_config.additive.mixture_scales = set_default_scale_prior(type);
    }

    if (auto dscale = fit.present<std::vector<double>>("--dscale"))
    {
        prior_config.dominant.mixture_scales
            = Eigen::Map<const Eigen::VectorXd>(
                dscale->data(), static_cast<Eigen::Index>(dscale->size()));
    }
    else
    {
        prior_config.dominant.mixture_scales = set_default_scale_prior(type);
    }
}

void setup_parallelization(
    int num_threads,
    const std::shared_ptr<spdlog::logger>& logger)
{
    if (num_threads > 0)
    {
        omp_set_num_threads(num_threads);
        Eigen::setNbThreads(num_threads);
        logger->info("Using {} threads for parallel computation", num_threads);
    }
}

void add_effect_to_model(
    gelex::BayesModel& model,
    gelex::GenotypeMap&& genotype_data,
    bool is_dominance)
{
    if (is_dominance)
    {
        model.add_dominance(std::move(genotype_data));
    }
    else
    {
        model.add_additive(std::move(genotype_data));
    }
}

void add_effect_to_model(
    gelex::BayesModel& model,
    gelex::GenotypeMatrix&& genotype_data,
    bool is_dominance)
{
    if (is_dominance)
    {
        model.add_dominance(std::move(genotype_data));
    }
    else
    {
        model.add_additive(std::move(genotype_data));
    }
}

int process_genotype_mmap(
    gelex::BayesModel& model,
    const std::shared_ptr<gelex::SampleManager>& sample_manager,
    const std::filesystem::path& bed_path,
    const std::string& out_prefix,
    int chunk_size,
    bool is_dominance,
    const std::shared_ptr<spdlog::logger>& logger)
{
    const std::string mmap_suffix = is_dominance ? ".dom" : ".add";
    const std::string mmap_path = out_prefix + mmap_suffix;

    auto genotype_pipe_result
        = gelex::GenotypePipe::create(bed_path, sample_manager, mmap_path);

    if (genotype_pipe_result)
    {
        auto process_result
            = is_dominance
                  ? genotype_pipe_result
                        ->process<gelex::DominantStandardizingProcessor>(
                            chunk_size)
                  : genotype_pipe_result
                        ->process<gelex::StandardizingProcessor>(chunk_size);
        VALIDATE_RESULT_OR_RETURN(process_result, logger);
        add_effect_to_model(
            model, std::move(process_result.value()), is_dominance);
    }
    else if (
        genotype_pipe_result.error().code == gelex::ErrorCode::OutputFileExists)
    {
        auto map_result = gelex::GenotypeMap::create(mmap_path + ".bmat");
        VALIDATE_RESULT_OR_RETURN(map_result, logger);
        add_effect_to_model(model, std::move(map_result.value()), is_dominance);
    }
    else
    {
        logger->error(genotype_pipe_result.error().message);
        return 1;
    }
    return 0;
}

int process_genotype_in_memory(
    gelex::BayesModel& model,
    const std::shared_ptr<gelex::SampleManager>& sample_manager,
    const std::filesystem::path& bed_path,
    int chunk_size,
    bool is_dominance,
    const std::shared_ptr<spdlog::logger>& logger)
{
    auto genotype_loader
        = gelex::GenotypeLoader::create(bed_path, sample_manager);
    VALIDATE_RESULT_OR_RETURN(genotype_loader, logger);

    auto process_result
        = is_dominance
              ? genotype_loader->process<gelex::DominantStandardizingProcessor>(
                    chunk_size)
              : genotype_loader->process<gelex::StandardizingProcessor>(
                    chunk_size);
    VALIDATE_RESULT_OR_RETURN(process_result, logger);

    add_effect_to_model(model, std::move(process_result.value()), is_dominance);
    return 0;
}

int process_genotype_effect(
    gelex::BayesModel& model,
    const std::shared_ptr<gelex::SampleManager>& sample_manager,
    const std::filesystem::path& bed_path,
    const std::string& out_prefix,
    int chunk_size,
    bool use_mmap,
    bool is_dominance,
    const std::shared_ptr<spdlog::logger>& logger)
{
    if (use_mmap)
    {
        return process_genotype_mmap(
            model,
            sample_manager,
            bed_path,
            out_prefix,
            chunk_size,
            is_dominance,
            logger);
    }

    return process_genotype_in_memory(
        model, sample_manager, bed_path, chunk_size, is_dominance, logger);
}

int configure_model_priors(
    gelex::BayesModel& model,
    gelex::BayesAlphabet type,
    argparse::ArgumentParser& fit,
    const std::shared_ptr<spdlog::logger>& logger)
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
    app::set_pi_prior(type, fit, prior_config);
    app::set_scale_prior(type, fit, prior_config);

    auto prior_result = (*prior_strategy)(model, prior_config);
    VALIDATE_RESULT_OR_RETURN(prior_result, logger);

    return 0;
}

int run_mcmc_analysis(
    gelex::BayesModel& model,
    gelex::BayesAlphabet type,
    const gelex::MCMCParams& mcmc_params,
    const std::filesystem::path& bim_path,
    std::string_view out_prefix,
    const std::shared_ptr<spdlog::logger>& logger)
{
    auto run_and_write = [&](auto trait_model)
    {
        gelex::MCMC mcmc(mcmc_params, trait_model);
        gelex::MCMCResult result = mcmc.run(model);
        gelex::MCMCResultWriter writer(result, bim_path);
        writer.save(out_prefix);
    };

    using bt = gelex::BayesAlphabet;
    switch (type)
    {
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
            run_and_write(gelex::BayesBd{});
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
            run_and_write(gelex::BayesRRd{});
            break;
        default:
            logger->error("Unsupported method: {}", type);
            return 1;
    }
    return 0;
}
}  // namespace app
