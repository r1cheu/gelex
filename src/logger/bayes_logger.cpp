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

void MCMCLogger::log_model_information(const BayesModel& model)
{
    logger_->info("");
    logger_->info("Model specification:");
    std::string genetic_terms;
    if (const auto* effect = model.additive(); effect != nullptr)
    {
        genetic_terms += " + add";
    }
    if (const auto* effect = model.dominant(); effect != nullptr)
    {
        genetic_terms += " + dom";
    }
    logger_->info("Priors:");

    if (const auto& effects = model.random(); !effects.empty())
    {
        for (const auto& effect : effects)
        {
            std::string name
                = effect.levels ? effect.levels.value()[0] : "test";
            logger_->info(
                "  {}: {}",
                name,
                sigma_prior("", effect.prior.nu, effect.prior.s2));
        }
    }

    if (const auto* effect = model.additive(); effect != nullptr)
    {
        logger_->info(
            "  add: {}", sigma_prior("", effect->prior.nu, effect->prior.s2));

        if (effect->pi.size() > 1)
        {
            std::string pi_str = "  π: [";
            for (Eigen::Index i = 0; i < effect->pi.size(); ++i)
            {
                if (i > 0)
                {
                    pi_str += ", ";
                }
                pi_str += fmt::format("{:.4f}", effect->pi(i));
            }
            pi_str += "]";
            logger_->info(pi_str);
        }
    }
    if (const auto* effect = model.dominant(); effect != nullptr)
    {
        logger_->info(
            "  dom_ratio: mean {:.3f}, var {:.3f}",
            effect->ratio_mean,
            effect->ratio_variance);
    }

    const auto& residual = model.residual();
    logger_->info(
        "  e: {}", sigma_prior("_e", residual.prior.nu, residual.prior.s2));

    logger_->info("");
    logger_->info("MCMC sampling started...");
}

void MCMCLogger::log_result(const MCMCResult& results, const BayesModel& model)
{
    logger_->info("");
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

    if (auto&& [effect, result]
        = std::make_pair(model.fixed(), results.fixed());
        effect != nullptr && result != nullptr)
    {
        if (effect->levels)
        {
            for (Eigen::Index i = 0;
                 i < static_cast<Eigen::Index>(effect->levels->size());
                 ++i)
            {
                log_summary(i, result->coeffs, effect->levels.value()[i]);
            }
        }
    }

    if (const auto* result = results.additive(); result != nullptr)
    {
        log_summary(0, result->variance, sigma_squared("_add"));
    }

    if (const auto* effect = model.additive();
        effect != nullptr && effect->pi.size() > 1)
    {
        if (const auto* pi_result = results.pi(); pi_result != nullptr)
        {
            for (Eigen::Index i = 0; i < pi_result->prop.size(); ++i)
            {
                log_summary(i, pi_result->prop, fmt::format("π[{}]", i));
            }
        }
    }
    if (const auto* result = results.dominant(); result != nullptr)
    {
        log_summary(0, result->variance, sigma_squared("_dom"));
    }
    log_summary(0, results.residual(), sigma_squared("_e"));
}

}  // namespace detail
}  // namespace gelex
