#include "gelex/estimator/bayes/logger.h"

#include <vector>

#include <barkeep.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
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
        fmt::join(model.random().names(), " + "),
        fmt::join(model.genetic().names(), " + ")));
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
            item("Random:  {}", fmt::join(model.random().names(), ", ")));
    }
    if (model.genetic())
    {
        logger_->info(
            item("Genetic:  {}", fmt::join(model.genetic().names(), ", ")));
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
}

std::shared_ptr<barkeep::CompositeDisplay> MCMCLogger::progress_bar(
    std::vector<size_t>& idxs,
    size_t total)
{
    std::vector<std::shared_ptr<bk::BaseDisplay>> displays;
    for (size_t i = 0; i < idxs.size(); ++i)
    {
        idxs[i] = 0;
        auto pb_config = bk::ProgressBarConfig<size_t>{
            .total = total,
            .message = fmt::format("Chain {}", i + 1),
            .speed = 0.1,
            .style = bk::ProgressBarStyle::Rich,
            .show = false};
        displays.push_back(bk::ProgressBar(&idxs[i], pb_config));
    }
    return bk::Composite(displays, "\n");
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
        header += format_names(model.random().names()) + "  ";
    }
    if (model.genetic())
    {
        header += format_names(model.genetic().names()) + "  ";
    }
    header += sigma_squared("_e") + "  h\u00B2  ETA";
    logger_->info(header);
}

void MCMCLogger::log_result(const MCMCResult& result)
{
    logger_->info(title(" Summary "));
    logger_->info(subtitle("Postier Results"));
    logger_->info(
        item("\u03BC: {}", with_std(result.mu.mean(0), result.mu.std(0))));
    logger_->info(item(
        "e: {}", with_std(result.residual.mean(0), result.residual.std(0))));
};
}  // namespace gelex
