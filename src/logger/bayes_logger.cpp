#include "../src/logger/bayes_logger.h"

#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>
#include "barkeep.h"

#include "../utils/formatter.h"
#include "gelex/logger.h"
#include "gelex/model/bayes/model.h"

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
    logger_->info(gelex::section("Model specification (Priors)..."));

    if (const auto& effects = model.random(); !effects.empty())
    {
        for (const auto& effect : effects)
        {
            std::string name
                = effect.levels ? effect.levels.value()[0] : "test";
            logger_->info(gelex::task("{}(rand)", name));
            logger_->info(
                gelex::subtask(
                    "Variance: Scaled Inv-χ²(ν={:.4f}, S²={:.4f}), init: "
                    "{:.4f}",
                    effect.prior.nu,
                    effect.prior.s2,
                    effect.init_variance));
        }
    }
    auto log_effect = [this](const auto* effect, bool dom)
    {
        if (!effect)
        {
            return;
        }
        std::string label = dom ? "Dominance" : "Additive";
        logger_->info(gelex::task("{} effect:", label));
        logger_->info(
            gelex::subtask(
                "Variance: Scaled Inv-χ²(ν={:.4f}, S²={:.4f}), init: {:.4f}",
                effect->marker_variance_prior.nu,
                effect->marker_variance_prior.s2,
                effect->init_marker_variance));

        if (effect->init_pi && effect->init_pi->size() > 1)
        {
            std::string pi_str = "[";
            for (Eigen::Index i = 0; i < effect->init_pi->size(); ++i)
            {
                if (i > 0)
                {
                    pi_str += ", ";
                }
                pi_str += fmt::format("{:.2f}", (*effect->init_pi)(i));
            }
            pi_str += "]";
            logger_->info(gelex::subtask("Mixture: {}", pi_str));
        }
    };

    log_effect(model.additive(), false);
    log_effect(model.dominant(), true);

    const auto& residual = model.residual();
    logger_->info(gelex::task("Residual:"));
    logger_->info(
        gelex::subtask(
            "Variance: Scaled Inv-χ²(ν={:.4f}, S²={:.4f}), init: {:.4f}",
            residual.prior.nu,
            residual.prior.s2,
            residual.init_variance));

    logger_->info("");
    logger_->info(gelex::section("MCMC Sampling..."));
}

void MCMCLogger::log_result(
    const MCMCResult& results,
    const BayesModel& model,
    double elapsed_time,
    Eigen::Index samples_collected)
{
    logger_->info("");
    logger_->info(
        fmt::format(
            fmt::emphasis::bold | fmt::fg(fmt::color::light_cyan),
            "── Posterior Summary "
            "─────────────────────────────────────────────────"));

    logger_->info("  Time elapsed: {:.2f}s", elapsed_time);
    logger_->info("  Samples collected per parameter: {}", samples_collected);
    logger_->info("");

    std::vector<std::string> header{
        "Parameter", "Mean", "SD", "5%", "95%", "n_eff", "r_hat"};

    auto format_header = [](const std::vector<std::string>& header)
    {
        std::string result = fmt::format("  {:<8} ", header[0]);
        for (size_t i = 1; i < header.size(); ++i)
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

    logger_->info(
        "  ───────────────────"
        "──────────────────────────────────────────────");

    auto log_summary = [this](
                           Eigen::Index i,
                           const PosteriorSummary& summary,
                           const std::string& name)
    {
        if (summary.rhat(i) > 1.1)
        {
            logger_->info(
                "  {:<8} {:>8.4f} {:>8.4f} {:>8.4f} {:>8.4f} "
                "{:>8.4f} {:>8.4f}*",
                name,
                summary.mean(i),
                summary.stddev(i),
                summary.hpdi_low(i),
                summary.hpdi_high(i),
                summary.ess(i),
                summary.rhat(i));
        }
        else
        {
            logger_->info(
                "  {:<8} {:>8.4f} {:>8.4f} {:>8.4f} {:>8.4f} "
                "{:>8.4f} {:>8.4f}",
                name,
                summary.mean(i),
                summary.stddev(i),
                summary.hpdi_low(i),
                summary.hpdi_high(i),
                summary.ess(i),
                summary.rhat(i));
        }
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
        log_summary(0, result->variance, "σ²_add");
        log_summary(0, result->heritability, "h²");
    }
    log_mixture(model.additive(), results.additive());

    if (const auto* result = results.dominant(); result != nullptr)
    {
        log_summary(0, result->variance, "σ²_dom");
        log_summary(0, result->heritability, "δ²");
    }
    log_mixture(model.dominant(), results.dominant());

    log_summary(0, results.residual(), "σ²_e");
    logger_->info(
        "  ───────────────────"
        "──────────────────────────────────────────────");
    logger_->info(
        "  * Values with high R-hat (>1.1) indicating poor convergence.");
    logger_->info("");
}

}  // namespace detail
}  // namespace gelex
