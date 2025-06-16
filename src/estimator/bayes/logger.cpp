#include "gelex/estimator/bayes/logger.h"

#include <vector>

#include <barkeep.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <armadillo>

#include "gelex/model/bayes/model.h"
#include "gelex/utils/formatter.h"
#include "gelex/utils/utils.h"

namespace gelex
{
namespace bk = barkeep;
/**
 * @brief Format the scale inverse chi-squared distribution
 * parameters into a human-readable string (e.g., "νS²χ⁻²(ν = nu, S² = s2)").
 *
 * @param nu Degrees of freedom.
 * @param s2 Scale parameter.
 * @return A formatted string representing the distribution parameters.
 */
std::string format_invchi(double nu, double s2)
{
    return fmt::format(
        "\u03BDS\u00B2\u03C7\u207B\u00B2(\u03BD = {}, S\u00B2 = {})", nu, s2);
}

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
    logger_->info(
        " \u25AA formula:  {}{}{} + e",
        model.formula(),
        fmt::join(model.random().names(), " +"),
        fmt::join(model.genetic().names(), " +"));

    logger_->info(" \u25AA Samples:  {:d}", model.n_individuals());
    logger_->info(" \u25AA Iters:  {:d}", params.n_iters);
    logger_->info(" \u25AA Burn-in:  {:d}", params.n_burnin);
    logger_->info(" \u25AA Thinning:  {:d}", params.n_thin);
    logger_->info(" \u25AA Chains:  {:d}", params.n_chains);
    logger_->info("");

    logger_->info(subtitle("Term Summary"));
    if (model.fixed())
    {
        logger_->info(
            " \u25AA Fixed:  {}", fmt::join(model.fixed()->names, ", "));
    };
    if (model.random())
    {
        logger_->info(
            " \u25AA Random:  {}", fmt::join(model.random().names(), ", "));
    }
    if (model.genetic())
    {
        logger_->info(
            " \u25AA Genetic:  {}", fmt::join(model.genetic().names(), ", "));
    }
    logger_->info("");

    logger_->info(subtitle("Prior Summary"));
    logger_->info(
        " \u25AA {} \u223C {}",
        format_sigma_squared(model.residual().name),
        format_invchi(model.residual().prior.nu, model.residual().prior.s2));
    logger_->info(title(" MCMC "));
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
            formatted.push_back(format_sigma_squared(name));
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
    header += format_sigma_squared("e") + "  h\u00B2  ETA";
    logger_->info(header);
}

void MCMCLogger::log_result(const MCMCResult& result)
{
    logger_->info(title(" Summary "));
    logger_->info(subtitle("Postier Results"));
    logger_->info(
        " \u25AA \u03BC: {}",
        format_value_with_std(result.mu.mean(0), result.mu.std(0)));
    logger_->info(
        " \u25AA {}: {}",
        format_sigma_squared("e"),
        format_value_with_std(result.residual.mean(0), result.residual.std(0)));
};
}  // namespace gelex
