#include "logger/bayes_logger.h"

#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <gelex/barkeep.h>

#include "gelex/logger.h"
#include "gelex/model/bayes/model.h"
#include "logger/logger_utils.h"

namespace gelex
{
namespace bk = barkeep;

namespace detail
{

MCMCLogger::MCMCLogger() : logger_{gelex::logging::get()} {}

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
    logger_->info("");
    logger_->info("Bayesian analysis started...");
    logger_->info("Model specification:");
    std::string genetic_terms;
    if (model.additive())
    {
        genetic_terms += " + " + model.additive()->name;
    }
    if (model.dominant())
    {
        genetic_terms += " + " + model.dominant()->name;
    }

    logger_->info(
        " Iterations: {:d}, Burn-in: {:d}, Thinning: {:d}, Chains: {:d}",
        params.n_iters,
        params.n_burnin,
        params.n_thin,
        params.n_chains);

    logger_->info("Priors:");

    if (model.random())
    {
        for (const auto& effect : model.random().effects())
        {
            logger_->info(
                "  {}: {}",
                effect.name,
                sigma_prior("", effect.prior.nu, effect.prior.s2));
        }
    }

    if (model.additive() && model.trait())
    {
        const auto& effect = *model.additive();
        auto prior_str = model.trait()->prior_info(
            effect.prior.nu, effect.prior.s2, effect.pi);
        for (size_t i{}; i < prior_str.size(); ++i)
        {
            if (i == 0)
            {
                logger_->info("  {}: {}", effect.name, prior_str[i]);
            }
            else
            {
                logger_->info("  {}", prior_str[i]);
            }
        }
    }
    if (model.dominant())
    {
        const auto& effect = *model.dominant();
        logger_->info(
            "  {}: mean {:.3f}, var {:.3f}",
            effect.name,
            effect.prior_mean,
            effect.prior_var);
    }

    logger_->info(
        "  e: {}",
        sigma_prior(
            "_e", model.residual().prior.nu, model.residual().prior.s2));

    logger_->info("");
    logger_->info("MCMC sampling started...");
}

void MCMCLogger::log_result(const MCMCResult& result, const BayesModel& model)
{
    logger_->info("MCMC results summary:");

    std::vector<std::string> header{
        "", "mean", "std", "5%", "95%", "n_eff", "r_hat"};

    auto format_header = [](const std::vector<std::string>& header)
    {
        std::string result;
        for (size_t i = 0; i < header.size(); ++i)
        {
            result += fmt::format("{:>8}", header[i]);
            if (i + 1 < header.size())
            {
                result += " ";
            }
        }
        return result;
    };
    logger_->info(format_header(header));

    auto log_summary = [this](
                           Eigen::Index i,
                           const PosteriorSummary& summary,
                           const std::string& name)
    {
        logger_->info(
            "{:>8} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} "
            "{:>8.1f} {:>8.2f}",
            name,
            summary.mean(i),
            summary.stddev(i),
            summary.hpdi_low(i),
            summary.hpdi_high(i),
            summary.ess(i),
            summary.rhat(i));
    };

    for (Eigen::Index i = 0;
         i < static_cast<Eigen::Index>(model.fixed().names.size());
         ++i)
    {
        log_summary(i, result.fixed(), model.fixed().names[i]);
    }

    if (model.additive())
    {
        const auto& effect = *model.additive();
        log_summary(
            0, result.additive()->variance, sigma_squared("_" + effect.name));
        if (model.trait() && model.trait()->estimate_pi())
        {
            // Handle pi estimation if the trait supports it
            // This would need access to pi results from the additive summary
        }
    }

    if (model.dominant())
    {
        const auto& effect = *model.dominant();
        log_summary(
            0, result.dominant()->variance, sigma_squared("_" + effect.name));
    }
    log_summary(0, result.residual(), sigma_squared("_e"));
}

}  // namespace detail
}  // namespace gelex
