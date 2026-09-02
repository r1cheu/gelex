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

#ifndef GELEX_BAYES_DETAIL_RECIPE_VALIDATION_H_
#define GELEX_BAYES_DETAIL_RECIPE_VALIDATION_H_

#include <string>
#include <vector>

#include "gelex/bayes/genetic/independent_topology.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/variance_budget.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::detail
{

// Gathers every domain violation of one recipe's inputs before reporting, so
// that a caller fixes all of them in one pass. Only throw_if_any() throws.
//
// Each issue is labelled with the modes it belongs to, because a leaf spec
// cannot know where it sits: the container holding it passes the scope in.
class RecipeIssues
{
   public:
    auto add(GeneticModeSet scope, std::string issue) -> void;
    auto throw_if_any() const -> void;

   private:
    std::vector<std::string> issues_;
};

auto check(RecipeIssues& issues, const SpikeSlab& spec, GeneticModeSet scope)
    -> void;
auto check(
    RecipeIssues& issues,
    const ScaledMixture& spec,
    GeneticModeSet scope) -> void;
auto check(
    RecipeIssues& issues,
    const JointSpikeSlab& spec,
    GeneticModeSet scope) -> void;
auto check(
    RecipeIssues& issues,
    const VarianceBudget& budget,
    GeneticModeSet modes) -> void;

constexpr auto check(
    RecipeIssues& /*issues*/,
    const NoParameters& /*spec*/,
    GeneticModeSet /*scope*/) noexcept -> void
{
}

// The topology already knows which mode each leaf stands for, so it narrows the
// scope instead of passing the recipe's own one down.
template <GeneticModeSet Modes, typename T>
auto check(
    RecipeIssues& issues,
    const IndependentTopology<Modes, T>& topology,
    GeneticModeSet /*scope*/) -> void
{
    for (const auto [mode, spec] : topology.each())
    {
        check(issues, spec, GeneticModeSet{mode});
    }
}

template <typename Parameters>
auto validate_recipe_inputs(
    const Parameters& parameters,
    const VarianceBudget& variance,
    GeneticModeSet modes) -> void
{
    auto issues = RecipeIssues{};

    check(issues, parameters, modes);
    check(issues, variance, modes);

    issues.throw_if_any();
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_DETAIL_RECIPE_VALIDATION_H_
