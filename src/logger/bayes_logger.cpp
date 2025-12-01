#include "logger/bayes_logger.h"

#include <vector>

#include <barkeep.h>
#include <fmt/format.h>
#include <fmt/ranges.h>

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
    auto log_effect = [this](const auto* effect, const char* label)
    {
        if (!effect)
        {
            return;
        }
        logger_->info(
            "  {}: {}, init: {:.6f}",
            label,
            sigma_prior(
                "",
                effect->marker_variance_prior.nu,
                effect->marker_variance_prior.s2),
            effect->init_marker_variance);

        if (effect->init_pi && effect->init_pi->size() > 1)
        {
            std::string pi_str = "  π: [";
            for (Eigen::Index i = 0; i < effect->init_pi->size(); ++i)
            {
                if (i > 0)
                {
                    pi_str += ", ";
                }
                pi_str += fmt::format("{:.4f}", (*effect->init_pi)(i));
            }
            pi_str += "]";
            logger_->info(pi_str);
        }
    };

    log_effect(model.additive(), "add");
    log_effect(model.dominant(), "dom");

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
            "{:>8} {:>8.4f} {:>8.4f} {:>8.4f} {:>8.4f} "
            "{:>8.4f} {:>8.4f}",
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

    auto log_mixture = [&](const auto* effect, const auto* result)
    {
        if (effect != nullptr && effect->init_pi && effect->init_pi->size() > 1)
        {
            if (result != nullptr && result->mixture_proportion.size() > 0)
            {
                for (Eigen::Index i = 0; i < result->mixture_proportion.size();
                     ++i)
                {
                    log_summary(
                        i, result->mixture_proportion, fmt::format("π[{}]", i));
                }
            }
        }
    };

    if (const auto* result = results.additive(); result != nullptr)
    {
        log_summary(0, result->variance, sigma_squared("_add"));
        log_summary(0, result->heritability, "h²");
    }
    log_mixture(model.additive(), results.additive());

    if (const auto* result = results.dominant(); result != nullptr)
    {
        log_summary(0, result->variance, sigma_squared("_dom"));
        log_summary(0, result->heritability, "δ²");
    }
    log_mixture(model.dominant(), results.dominant());

    log_summary(0, results.residual(), sigma_squared("_e"));
}

}  // namespace detail
}  // namespace gelex
