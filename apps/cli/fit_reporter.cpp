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
#include <string_view>
#include <type_traits>
#include <variant>

#include <fmt/format.h>
#include <Eigen/Core>

#include "cli/report_printer.h"
#include "gelex/algo/infer/posterior_summary.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/parameter/values.h"
#include "gelex/bayes/prior.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/logging/formatter.h"
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
    const bayes::ScaledInvChiSqPrior& prior,
    double init_variance) -> void
{
    cli::printer().line(
        "    Variance: Scaled Inv-χ²(ν={:.4f}, S²={:.4f}), init: {:.4f}",
        prior.degrees_of_freedom(),
        prior.scale(),
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

auto FitReporter::print_random_prior(const bayes::RandomPrior& prior) -> void
{
    cli::printer().line("   Random effect:");
    print_variance_prior(prior.prior(), prior.initial_value());
}

auto FitReporter::print_genetic_prior(const bayes::GeneticPrior& prior) -> void
{
    auto format_vec = [](const auto& p)
    {
        auto formatted = std::span(p.data(), p.size())
                         | std::views::transform(
                             [](double v) { return fmt::format("{:.3f}", v); });
        return fmt::join(formatted, ", ");
    };

    auto print_proportion = [&format_vec](const auto& proportion)
    {
        using Proportion = std::remove_cvref_t<decltype(proportion)>;
        if constexpr (std::is_same_v<Proportion, bayes::FixedProportion>)
        {
            cli::printer().line(
                "    Proportion: [{}]", format_vec(proportion.value()));
        }
        else
        {
            cli::printer().line(
                "    Proportion: [{}]",
                format_vec(proportion.parameter().initial_value()));
        }
    };

    std::visit(
        [&](const auto& genetic_group)
        {
            using Group = std::remove_cvref_t<decltype(genetic_group)>;
            if constexpr (std::is_same_v<Group, bayes::SingleGeneticPrior>)
            {
                std::visit(
                    [&](const auto& genetic)
                    {
                        cli::printer().line("   {} effect:", genetic.mode());
                        const auto& parameter = genetic.variance().parameter();
                        print_variance_prior(
                            parameter.prior(), parameter.initial_value());
                        if constexpr (requires { genetic.proportion(); })
                        {
                            print_proportion(genetic.proportion());
                        }
                        if constexpr (requires { genetic.multiplier(); })
                        {
                            cli::printer().line(
                                "    Multiplier: [{}]",
                                format_vec(genetic.multiplier()));
                        }
                    },
                    genetic_group);
            }
            else
            {
                cli::printer().line("   A, D effect:");
                std::visit(
                    [&](const auto& genetic)
                    {
                        for (const auto mode : {GeneticMode::A, GeneticMode::D})
                        {
                            const auto& parameter
                                = genetic.variance(mode).parameter();
                            print_variance_prior(
                                parameter.prior(), parameter.initial_value());
                        }
                        print_proportion(genetic.proportion());
                    },
                    genetic_group);
            }
        },
        prior);
}

auto FitReporter::print_residual_prior(const bayes::ResidualPrior& prior)
    -> void
{
    cli::printer().line("   Residual:");
    print_variance_prior(prior.prior(), prior.initial_value());
}

}  // namespace gelex::cli
