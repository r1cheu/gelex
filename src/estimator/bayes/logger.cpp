#include "gelex/estimator/bayes/logger.h"

#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <gelex/barkeep.h>
#include <armadillo>

#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/policy.h"
#include "gelex/utils/formatter.h"
#include "gelex/utils/utils.h"

namespace gelex
{
namespace bk = barkeep;

MCMCLogger::MCMCLogger() : logger_{Logger::logger()} {}

void MCMCLogger::set_verbose(bool verbose)
{
    if (!verbose)
    {
        logger_->set_level(spdlog::level::warn);
    }
}

void MCMCLogger::log_model_information(
    const BayesModel& model,
    MCMCParams params)
{
    logger_->info(title(" BAYESALPHABET ANALYSIS "));
    logger_->info(subtitle("Model Specification"));
    logger_->info(item(
        "formula:  {}{}{} + e",
        model.formula(),
        fmt::join(model.random().keys(), " + "),
        fmt::join(model.genetic().keys(), " + ")));
    logger_->info(item("Samples:  {:d}", model.n_individuals()));
    logger_->info(item("Iters:  {:d}", params.n_iters));
    logger_->info(item("Burn-in:  {:d}", params.n_burnin));
    logger_->info(item("Thinning:  {:d}", params.n_thin));
    logger_->info(item("Chains:  {:d}", params.n_chains));
    logger_->info("");

    logger_->info(subtitle("Term Summary"));
    logger_->info(item("Fixed:  {}", fmt::join(model.fixed().names, ", ")));
    if (model.random())
    {
        logger_->info(
            item("Random:  {}", fmt::join(model.random().keys(), ", ")));
    }
    if (model.genetic())
    {
        logger_->info(
            item("Genetic:  {}", fmt::join(model.genetic().keys(), ", ")));
    }
    logger_->info("");

    logger_->info(subtitle("Prior Summary"));

    if (model.random())
    {
        for (const auto& effect : model.random().effects())
        {
            logger_->info(item(
                "{}: {}",
                effect.name,
                sigma_prior("", effect.prior.nu, effect.prior.s2)));
        }
    }

    if (model.genetic())
    {
        for (const auto& effect : model.genetic().effects())
        {
            auto prior_str = bayes_trait_prior_str[to_index(effect.type)](
                effect.prior.nu, effect.prior.s2, effect.pi);
            for (size_t i{}; i < prior_str.size(); ++i)
            {
                if (i == 0)
                {
                    logger_->info(item("{}: {}", effect.name, prior_str[i]));
                }
                else
                {
                    logger_->info("{}", prior_str[i]);
                }
            }
        }
    }

    logger_->info(item(
        "e: {}",
        sigma_prior(
            "_e", model.residual().prior.nu, model.residual().prior.s2)));

    logger_->info(title(" MCMC "));
}

void MCMCLogger::log_iter_header(const BayesModel& model)
{
    // Helper function to format names with sigma squared
    auto format_names = [](const std::vector<std::string>& names)
    {
        std::vector<std::string> formatted;
        formatted.reserve(names.size());
        for (const auto& name : names)
        {
            formatted.push_back(sigma_squared("_" + name));
        }
        return fmt::to_string(fmt::join(formatted, "  "));
    };

    std::string header = "Iter.  \u03BC  ";

    if (model.random())
    {
        header += format_names(model.random().keys()) + "  ";
    }
    if (model.genetic())
    {
        header += format_names(model.genetic().keys()) + "  ";
    }
    header += sigma_squared("_e") + "  h\u00B2  ETA";
    logger_->info(header);
}

void MCMCLogger::log_result(const MCMCResult& result, const BayesModel& model)
{
    logger_->info(title(" Summary "));

    std::vector<std::string> header{
        "", "mean", "std", "median", "5%", "95%", "n_eff", "r_hat"};

    auto format_header = [](const std::vector<std::string>& header)
    {
        std::string result;
        for (size_t i = 0; i < header.size(); ++i)
        {
            result += fmt::format("{:>9}", header[i]);
            if (i + 1 < header.size())
            {
                result += "  ";
            }
        }
        return result;
    };
    logger_->info(format_header(header));

    auto log_summary =
        [this](
            size_t i, const PosteriorSummary& summary, const std::string& name)
    {
        logger_->info(
            "{:>9}  {:>9.2f}  {:>9.2f}  {:>9.2f}  {:>9.2f}  {:>9.2f}  "
            "{:>9.2f}  {:>9.2f}",
            name,
            summary.mean.at(i),
            summary.stddev.at(i),
            summary.median.at(i),
            summary.hpdi_low.at(i),
            summary.hpdi_high.at(i),
            summary.ess.at(i),
            summary.rhat.at(i));
    };

    for (size_t i = 0; i < model.fixed().levels.size(); ++i)
    {
        log_summary(i, result.fixed, model.fixed().levels[i]);
    }

    for (size_t i = 0; i < model.genetic().size(); ++i)
    {
        log_summary(
            i,
            result.genetic[i].genetic_var,
            sigma_squared("_" + model.genetic()[i].name));
        log_summary(
            i,
            result.genetic[i].heritability,
            h2("_" + model.genetic()[i].name));

        if (bayes_trait_estimate_pi[to_index(model.genetic()[i].type)])
        {
            for (size_t j = 0; j < result.genetic[i].pi.mean.n_elem; ++j)
            {
                log_summary(
                    j,
                    result.genetic[i].pi,
                    fmt::format("π{}_{}", j, model.genetic()[i].name));
            }
        }
    }
    log_summary(0, result.residual, sigma_squared("_e"));
};
}  // namespace gelex
