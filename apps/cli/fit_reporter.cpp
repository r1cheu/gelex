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

#include "cli/report_printer.h"
#include "gelex/algo/infer/posterior_summary.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/logging/formatter.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/capabilities.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::cli
{

auto FitReporter::on_event(const FitPriorSetEvent& event) const -> void
{
    cli::printer().block(gelex::section("[Prior Configuration]"));

    if (event.prior == nullptr)
    {
        return;
    }
    const auto& prior = *event.prior;
    print_random_prior(prior.random());
    for (const auto& genetic : prior.genetics())
    {
        print_genetic_prior(genetic);
    }
    print_residual_prior(prior.residual());
}

auto FitReporter::on_event(const FitResultsSavedEvent& event) const -> void
{
    cli::printer().block(
        gelex::success(
            "Results saved to '{}' (.param, .summary, .snp.eff, .log)",
            event.out_prefix));
}

auto FitReporter::print_variance_prior(
    const stats::detail::ScaledInvChiSqParams& prior,
    double init_variance) -> void
{
    cli::printer().line(
        "    Variance: Scaled Inv-χ²(ν={:.4f}, S²={:.4f}), init: {:.4f}",
        prior.nu,
        prior.s2,
        init_variance);
}

auto FitReporter::print_summary_row(
    std::string_view name,
    const PosteriorSummary& summary,
    Eigen::Index index) -> void
{
    cli::printer().line(
        "  {:<8} {:>10.6f} {:>10.6f}",
        name,
        summary.mean(index),
        summary.stddev(index));
}

auto FitReporter::print_random_prior(const bayes::VarianceParameter& parameter)
    -> void
{
    cli::printer().line("   Random effect:");
    print_variance_prior(
        stats::detail::ScaledInvChiSqParams{
            parameter.prior().degrees_of_freedom(), parameter.prior().scale()},
        parameter.initial_value());
}

auto FitReporter::print_genetic_prior(const bayes::GeneticPrior& prior) -> void
{
    cli::printer().line("   {} effect:", fmt::join(prior.modes(), ", "));

    auto format_vec = [](const auto& p)
    {
        auto formatted = std::span(p.data(), p.size())
                         | std::views::transform(
                             [](double v) { return fmt::format("{:.3f}", v); });
        return fmt::join(formatted, ", ");
    };

    if (const auto* variance = prior.query<bayes::MarkerVarianceCap>())
    {
        for (const auto& marker_variance : variance->variance())
        {
            const auto& parameter = marker_variance.parameter();
            print_variance_prior(
                stats::detail::ScaledInvChiSqParams{
                    parameter.prior().degrees_of_freedom(),
                    parameter.prior().scale()},
                parameter.initial_value());
        }
    }
    if (const auto* proportion = prior.query<bayes::MixtureProportionCap>())
    {
        for (const auto& mixture_proportion : proportion->proportion())
        {
            cli::printer().line(
                "    Proportion: [{}]",
                format_vec(mixture_proportion.parameter().initial_value()));
        }
    }
    if (const auto* multiplier = prior.query<bayes::MultiplierCap>())
    {
        for (const auto& value : multiplier->multiplier())
        {
            cli::printer().line("    Multiplier: [{}]", format_vec(value));
        }
    }
}

auto FitReporter::print_residual_prior(
    const bayes::VarianceParameter& parameter) -> void
{
    cli::printer().line("   Residual:");
    print_variance_prior(
        stats::detail::ScaledInvChiSqParams{
            parameter.prior().degrees_of_freedom(), parameter.prior().scale()},
        parameter.initial_value());
}

}  // namespace gelex::cli
