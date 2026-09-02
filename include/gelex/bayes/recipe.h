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

#ifndef GELEX_BAYES_RECIPE_H_
#define GELEX_BAYES_RECIPE_H_

#include <concepts>
#include <utility>

#include "gelex/bayes/detail/recipe_validation.h"
#include "gelex/bayes/genetic/independent_topology.h"
#include "gelex/bayes/semantic_method.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/variance_budget.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{

namespace detail
{

template <typename Method, GeneticModeSet Modes>
struct MethodParameters;

template <Variance Kind, GeneticModeSet Modes>
struct MethodParameters<GaussianMethod<Kind>, Modes>
{
    using type = NoParameters;
};

template <Variance Kind, GeneticModeSet Modes>
struct MethodParameters<SpikeSlabMethod<Kind>, Modes>
{
    using type = IndependentTopology<Modes, SpikeSlab>;
};

template <GeneticModeSet Modes>
struct MethodParameters<ScaledMixtureMethod, Modes>
{
    using type = IndependentTopology<Modes, ScaledMixture>;
};

template <>
struct MethodParameters<JointSpikeSlabMethod, GeneticMode::A | GeneticMode::D>
{
    using type = JointSpikeSlab;
};

template <typename Method, GeneticModeSet Modes>
using method_parameters_t = typename MethodParameters<Method, Modes>::type;

template <typename SemanticMethod, GeneticModeSet Modes>
concept SupportedSemanticMethod = Modes.size() > 0 && requires {
    typename method_parameters_t<SemanticMethod, Modes>;
};

}  // namespace detail

template <typename SemanticMethod, GeneticModeSet Modes>
    requires detail::SupportedSemanticMethod<SemanticMethod, Modes>
class BayesRecipe
{
    using MethodParameterT = detail::method_parameters_t<SemanticMethod, Modes>;

   public:
    static constexpr GeneticModeSet modes = Modes;
    using method_type = SemanticMethod;

    BayesRecipe(MethodParameterT parameters, VarianceBudget variance)
        : parameters_(std::move(parameters)), variance_(variance)
    {
        validate();
    }

    explicit BayesRecipe(VarianceBudget variance)
        requires std::same_as<MethodParameterT, NoParameters>
        : variance_(variance)
    {
        validate();
    }

    [[nodiscard]] static auto defaults() -> BayesRecipe
    {
        return BayesRecipe{
            MethodParameterT{}, VarianceBudget{default_shares(Modes)}};
    }

    [[nodiscard]] auto parameters() const noexcept -> const MethodParameterT&
    {
        return parameters_;
    }

    [[nodiscard]] auto variance() const noexcept -> const VarianceBudget&
    {
        return variance_;
    }

   private:
    auto validate() const -> void
    {
        detail::validate_recipe_inputs(parameters_, variance_, Modes);
    }

    [[no_unique_address]] MethodParameterT parameters_;
    VarianceBudget variance_;
};

}  // namespace gelex

#endif  // GELEX_BAYES_RECIPE_H_
