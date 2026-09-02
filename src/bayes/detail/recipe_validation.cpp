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

#include "gelex/bayes/detail/recipe_validation.h"

#include <cmath>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <utility>

#include "gelex/bayes/spec.h"
#include "gelex/bayes/variance_budget.h"
#include "gelex/exception.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::detail
{

namespace
{

constexpr double SIMPLEX_TOLERANCE = 1e-9;

[[nodiscard]] auto is_finite_positive(double value) noexcept -> bool
{
    return std::isfinite(value) && value > 0.0;
}

auto check_simplex(
    RecipeIssues& issues,
    GeneticModeSet scope,
    std::span<const double> weights) -> void
{
    auto total = 0.0;
    for (const auto [index, weight] : weights | std::views::enumerate)
    {
        if (!is_finite_positive(weight))
        {
            issues.add(
                scope,
                fmt::format(
                    "mixture weights[{}] must be a finite positive "
                    "probability, got {}",
                    index,
                    weight));
        }
        total += weight;
    }

    if (std::abs(total - 1.0) > SIMPLEX_TOLERANCE)
    {
        issues.add(
            scope, fmt::format("mixture weights must sum to 1, got {}", total));
    }
}

auto check_open_unit(
    RecipeIssues& issues,
    GeneticModeSet scope,
    double value,
    std::string_view field) -> void
{
    if (!is_finite_positive(value) || value >= 1.0)
    {
        issues.add(
            scope,
            fmt::format(
                "{} must lie in the open interval (0, 1), got {}",
                field,
                value));
    }
}

}  // namespace

auto RecipeIssues::add(GeneticModeSet scope, std::string issue) -> void
{
    issues_.push_back(
        scope == GeneticModeSet{} ? std::move(issue)
                                  : fmt::format("{} {}", scope, issue));
}

auto RecipeIssues::throw_if_any() const -> void
{
    if (issues_.empty())
    {
        return;
    }

    throw GelexException(
        fmt::format(
            "invalid Bayes recipe input:\n  {}", fmt::join(issues_, "\n  ")));
}

auto check(RecipeIssues& issues, const SpikeSlab& spec, GeneticModeSet scope)
    -> void
{
    check_open_unit(
        issues, scope, spec.probability.initial, "inclusion probability");
}

auto check(
    RecipeIssues& issues,
    const ScaledMixture& spec,
    GeneticModeSet scope) -> void
{
    check_simplex(issues, scope, spec.probabilities.initial);

    // The first component is the null scale, which zeroes the effect it
    // weights. The array has a fixed extent, so there is always one.
    if (spec.scales.front() != 0.0)
    {
        issues.add(
            scope,
            fmt::format(
                "variance multipliers[0] must be zero, got {}",
                spec.scales.front()));
    }

    for (const auto [index, scale] :
         std::span{spec.scales}.subspan(1) | std::views::enumerate)
    {
        if (!is_finite_positive(scale))
        {
            issues.add(
                scope,
                fmt::format(
                    "variance multipliers[{}] must be a finite positive "
                    "multiplier, got {}",
                    index + 1,
                    scale));
        }
    }
}

auto check(
    RecipeIssues& issues,
    const JointSpikeSlab& spec,
    GeneticModeSet scope) -> void
{
    check_simplex(issues, scope, spec.probabilities.initial);
    check_open_unit(
        issues,
        scope,
        spec.positive_probability.initial,
        "dominance positive-sign probability");
}

auto check(
    RecipeIssues& issues,
    const VarianceBudget& budget,
    GeneticModeSet modes) -> void
{
    auto allocated = budget.random();

    if (!std::isfinite(budget.random()) || budget.random() < 0.0)
    {
        issues.add(
            {},
            fmt::format(
                "random variance share must be finite and non-negative, got {}",
                budget.random()));
    }

    for (const auto mode : ALL_GENETIC_MODES)
    {
        const auto share = budget.share(mode);
        allocated += share;

        const auto scope = GeneticModeSet{mode};

        if (modes.contains(mode))
        {
            if (!is_finite_positive(share))
            {
                issues.add(
                    scope,
                    fmt::format(
                        "variance share must be finite and positive, got {}",
                        share));
            }
        }
        else if (share != 0.0)
        {
            issues.add(
                scope,
                fmt::format(
                    "variance share must be zero when the mode is absent, "
                    "got {}",
                    share));
        }
    }

    if (!(allocated < 1.0))
    {
        issues.add(
            {},
            fmt::format(
                "variance shares must sum to less than 1, got {}", allocated));
    }
}

}  // namespace gelex::detail
