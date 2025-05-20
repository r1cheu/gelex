#include "gelex/estimator/mcmc_logger.h"
#include <fmt/ranges.h>
#include <vector>

#include "armadillo"
#include "gelex/model/bayes.h"
#include "gelex/utils.h"

namespace gelex
{

MCMCLogger::MCMCLogger() : logger_{Logger::logger()} {}

void MCMCLogger::set_verbose(bool verbose)
{
    if (!verbose)
    {
        logger_->set_level(spdlog::level::warn);
    }
}

void MCMCLogger::log_model_information(const Bayes& model, MCMCParams params)
{
    logger_->info(
        "─────────────────────── BAYESALPHABET ANALYSIS "
        "───────────────────────");
    logger_->info("");
    logger_->info(fmt::format("{}", wine_red("[Model Specification]")));
    logger_->info(
        " \u25AA Model:  {}{}{} + e",
        model.formula(),
        fmt::join(model.random().names(), " +"),
        fmt::join(model.genetic().names(), " +"));

    logger_->info(" \u25AA Samples:  {:d}", model.n_individuals());
    logger_->info(" \u25AA Iters. {:d}", params.iter);
    logger_->info(" \u25AA Burn-in: {:d}", params.n_burnin);
    logger_->info(" \u25AA Thinning: {:d}", params.n_thin);
    logger_->info("");

    logger_->info(fmt::format("{}", wine_red("[Term Summary]")));
    logger_->info(" \u25AA Fixed:  {}", fmt::join(model.fixed().names, ", "));
    if (model.random().has_effects())
    {
        logger_->info(
            " \u25AA Random:  {}", fmt::join(model.random().names(), ", "));
    }
    if (model.genetic().has_effects())
    {
        logger_->info(
            " \u25AA Genetic:  {}", fmt::join(model.genetic().names(), ", "));
    }
    logger_->info("");

    logger_->info(fmt::format("{}", wine_red("[Prior Summary]")));
    logger_->info(
        " \u25AA \u03C3_{}\u00B2 \u223C \u03BDS\u00B2\u03C7\u207B\u00B2(\u03BD "
        "= {}, S\u00B2 = {})",
        model.residual().name,
        model.residual().prior.nu,
        model.residual().prior.s2);
    logger_->info("");

    logger_->info(
        "───────────────────────── MCMC BURNIN "
        "────────────────────────");
}

void MCMCLogger::log_iteration(
    size_t iter,
    const Bayes& model,
    std::string_view duartion)
{
    double sum_var = 0.0;
    for (const auto& eff : model.random())
    {
        sum_var += arma::as_scalar(eff.sigma);
    }

    std::vector<double> genetic_vars;
    for (const auto& eff : model.genetic())
    {
        double sigma = arma::var(eff.u);
        sum_var += sigma;
        genetic_vars.push_back(sigma);
    }

    sum_var += model.residual().value;

    std::vector<double> h2;
    h2.reserve(genetic_vars.size());
    for (double& sigma : genetic_vars)
    {
        h2.push_back(sigma / sum_var);
    }

    logger_->info(
        "{:d} {:.3f} {:.3f} {:.3f} ({})",
        iter,
        fmt::join(genetic_vars, ", "),
        model.residual().value,
        fmt::join(h2, ", "),
        duartion);
}

void MCMCLogger::log_burnin_finished()
{
    logger_->info(
        "───────────────────────── MCMC BURNIN "
        "FINISHED ────────────────────────");
}
}  // namespace gelex
