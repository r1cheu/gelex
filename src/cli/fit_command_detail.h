#ifndef GELEX_CLI_FIT_COMMAND_DETAIL_H_
#define GELEX_CLI_FIT_COMMAND_DETAIL_H_
#include <filesystem>

#include <Eigen/Core>

#include "gelex/estimator/bayes/params.h"
#include "gelex/model/effects.h"

#define VALIDATE_RESULT_OR_RETURN(result, logger) \
    if (!result)                                  \
    {                                             \
        logger->error(result.error().message);    \
        return 1;                                 \
    }

namespace spdlog
{
class logger;
}

namespace gelex
{
class SampleManager;
struct PriorConfig;
class BayesModel;
class GenotypeMap;
class GenotypeMatrix;
struct MCMCParams;

}  // namespace gelex

namespace argparse
{
class ArgumentParser;
}

namespace app
{
auto has_dominance(gelex::BayesAlphabet type) -> bool;

auto set_default_pi_prior(gelex::BayesAlphabet type) -> Eigen::VectorXd;

auto set_pi_prior(
    gelex::BayesAlphabet type,
    argparse::ArgumentParser& fit,
    gelex::PriorConfig& prior_config) -> void;

auto set_default_scale_prior(gelex::BayesAlphabet type) -> Eigen::VectorXd;

auto set_scale_prior(
    gelex::BayesAlphabet type,
    argparse::ArgumentParser& fit,
    gelex::PriorConfig& prior_config) -> void;

auto has_dominance_support(gelex::BayesAlphabet type) -> bool;

int process_genotype_effect(
    gelex::BayesModel& model,
    const std::shared_ptr<gelex::SampleManager>& sample_manager,
    const std::filesystem::path& bed_path,
    const std::string& out_prefix,
    int chunk_size,
    bool use_mmap,
    bool is_dominance,
    const std::shared_ptr<spdlog::logger>& logger);

int configure_model_priors(
    gelex::BayesModel& model,
    gelex::BayesAlphabet type,
    argparse::ArgumentParser& fit,
    const std::shared_ptr<spdlog::logger>& logger);

int run_mcmc_analysis(
    gelex::BayesModel& model,
    gelex::BayesAlphabet type,
    const gelex::MCMCParams& mcmc_params,
    const std::filesystem::path& bim_path,
    std::string_view out_prefix,
    const std::shared_ptr<spdlog::logger>& logger);
}  // namespace app

#endif  // GELEX_CLI_FIT_COMMAND_DETAIL_H_
