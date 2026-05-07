/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "fit_reporter.h"

#include <ranges>
#include <span>

#include <fmt/format.h>

#include "gelex/algo/infer/posterior_summary.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/logging/formatter.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/method.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::cli
{

FitReporter::FitReporter() : logger_(gelex::logging::get()) {}

auto FitReporter::on_event(const FitMethodSetEvent& event) const -> void
{
    logger_->info(gelex::section("[Prior Configuration]"));

    if (event.method == nullptr)
    {
        return;
    }
    const auto& method = *event.method;
    for (const auto& spec : method.randoms)
    {
        print_random_prior(spec);
    }
    for (const auto& prior : method.genetics)
    {
        std::visit(
            [&](const auto& spec)
            {
                using T = std::decay_t<decltype(spec)>;
                if constexpr (std::is_same_v<T, bayes::GeneticSpec>)
                {
                    print_genetic_prior(prior, spec.mode);
                }
                else
                {
                    print_genetic_prior(prior, spec.additive.mode);
                    print_genetic_prior(prior, spec.dominance.mode);
                }
            },
            prior.spec);
    }
    print_residual_prior(method.residual);
}

auto FitReporter::on_event(const FitResultsSavedEvent& event) const -> void
{
    logger_->info(
        gelex::success(
            "Results saved to '{}' (.param, .summary, .snp.eff, .log)",
            event.out_prefix));
}

auto FitReporter::print_variance_prior(
    const stats::detail::ScaledInvChiSqParams& prior,
    double init_variance) const -> void
{
    logger_->info(
        "    Variance: Scaled Inv-χ²(ν={:.4f}, S²={:.4f}), init: {:.4f}",
        prior.nu,
        prior.s2,
        init_variance);
}

auto FitReporter::print_summary_row(
    std::string_view name,
    const PosteriorSummary& summary,
    Eigen::Index index) const -> void
{
    logger_->info(
        "  {:<8} {:>10.6f} {:>10.6f}",
        name,
        summary.mean(index),
        summary.stddev(index));
}

auto FitReporter::print_random_prior(const bayes::VarianceSpec& spec) const
    -> void
{
    logger_->info("   Random effect:");
    print_variance_prior(
        stats::detail::ScaledInvChiSqParams{spec.prior.nu, spec.prior.s2},
        spec.init);
}

auto FitReporter::print_genetic_prior(
    const bayes::GeneticPrior& prior,
    GeneticMode mode) const -> void
{
    logger_->info("   {} effect:", mode);

    auto format_vec = [](const auto& p)
    {
        auto formatted = std::span(p.data(), p.size())
                         | std::views::transform(
                             [](double v) { return fmt::format("{:.3f}", v); });
        return fmt::join(formatted, ", ");
    };

    const auto& spec = std::visit(
        [&](const auto& s) -> const bayes::GeneticSpec&
        {
            using T = std::decay_t<decltype(s)>;
            if constexpr (std::is_same_v<T, bayes::GeneticSpec>)
            {
                return s;
            }
            else
            {
                return (s.additive.mode == mode) ? s.additive : s.dominance;
            }
        },
        prior.spec);

    print_variance_prior(
        stats::detail::ScaledInvChiSqParams{
            spec.variance.prior.nu, spec.variance.prior.s2},
        spec.variance.init);

    if (prior.mixture)
    {
        logger_->info(
            "    Proportion: [{}]",
            format_vec(prior.mixture->proportions.init));
        if (const auto* sm
            = std::get_if<bayes::ScaledMixture>(&prior.mixture->strategy))
        {
            logger_->info("    Multiplier: [{}]", format_vec(sm->multiplier));
        }
    }
}

auto FitReporter::print_residual_prior(const bayes::VarianceSpec& spec) const
    -> void
{
    logger_->info("   Residual:");
    print_variance_prior(
        stats::detail::ScaledInvChiSqParams{spec.prior.nu, spec.prior.s2},
        spec.init);
}

}  // namespace gelex::cli
