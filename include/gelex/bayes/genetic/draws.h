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

#ifndef GELEX_BAYES_GENETIC_DRAWS_H_
#define GELEX_BAYES_GENETIC_DRAWS_H_

#include <cstddef>
#include <type_traits>
#include <utility>

#include "gelex/bayes/basic_draw.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/genetic_mode.h"
#include "gelex/infra/var.h"

namespace gelex
{

namespace detail
{

template <VarianceLayout Kind>
struct MarkerVarianceDrawType;

template <>
struct MarkerVarianceDrawType<VarianceLayout::Pooled>
{
    using type = ScalarDraw;
};

template <>
struct MarkerVarianceDrawType<VarianceLayout::Unpooled>
{
    using type = VectorDraw;
};

template <VarianceLayout Kind>
using marker_variance_draw_t = typename MarkerVarianceDrawType<Kind>::type;

template <UpdatePolicy Policy, typename Draw>
using policy_draw_t
    = std::conditional_t<Policy == UpdatePolicy::Fixed, EmptyDraw, Draw>;

}  // namespace detail

template <VarianceLayout Kind>
struct GaussianDraws
{
    detail::marker_variance_draw_t<Kind> variance;

    auto append(const GaussianState<Kind>& state) -> void
    {
        variance.append(state.variance);
    }
};

template <HalfNormalAsymmetry Axis>
struct HalfNormalDraws;

template <>
struct HalfNormalDraws<HalfNormalAsymmetry::Count>
{
    ScalarDraw variance;
    CategoryDraw<3> assignment;
    ScalarDraw positive_probability;

    auto append(const HalfNormalState<HalfNormalAsymmetry::Count>& state)
        -> void
    {
        variance.append(state.variance);
        assignment.append(state.assignment);
        positive_probability.append(state.positive_probability);
    }
};

template <>
struct HalfNormalDraws<HalfNormalAsymmetry::Magnitude>
{
    VectorDraw variances;
    CategoryDraw<3> assignment;

    auto append(const HalfNormalState<HalfNormalAsymmetry::Magnitude>& state)
        -> void
    {
        variances.append(state.variances);
        assignment.append(state.assignment);
    }
};

template <VarianceLayout Kind, UpdatePolicy ProbabilityUpdate>
struct SpikeSlabDraws
{
    detail::marker_variance_draw_t<Kind> variance;
    CategoryDraw<2> assignment;
    detail::policy_draw_t<ProbabilityUpdate, ScalarDraw> probability;

    auto append(const SpikeSlabState<Kind>& state) -> void
    {
        variance.append(state.variance);
        assignment.append(state.assignment);
        probability.append(state.probability);
    }
};

template <std::size_t ClassCount, UpdatePolicy ProbabilitiesUpdate>
struct ScaledMixtureDraws
{
    ScalarDraw variance;
    CategoryDraw<ClassCount> assignment;
    detail::policy_draw_t<ProbabilitiesUpdate, VectorDraw> probabilities;
    VectorDraw component_explained_variance;

    auto append(const ScaledMixtureState<ClassCount>& state) -> void
    {
        variance.append(state.variance);
        assignment.append(state.assignment);
        probabilities.append(state.probabilities);
        component_explained_variance.append(
            matvar<0>(state.fitted_values, VarNormType::Population));
    }
};

template <std::size_t ClassCount, UpdatePolicy ProbabilitiesUpdate>
struct JointSpikeSlabDraws
{
    CategoryDraw<ClassCount> assignment;
    detail::policy_draw_t<ProbabilitiesUpdate, VectorDraw> probabilities;
    VectorDraw component_explained_variance;

    auto append(const JointSpikeSlabState<ClassCount>& state) -> void
    {
        assignment.append(state.assignment);
        probabilities.append(state.probabilities);
        component_explained_variance.append(
            matvar<0>(state.fitted_values, VarNormType::Population));
    }
};

template <typename FamilyDraws>
struct GeneticModeDraws
{
    VectorDraw coefficients;
    FamilyDraws family_draws;

    template <typename FamilyState>
    auto append(const GeneticModeState<FamilyState>& state) -> void
    {
        coefficients.append(state.coefficients);
        family_draws.append(state.family_state);
    }
};

template <typename GeneticState, typename ModeDraws>
class IndependentDraws
{
   public:
    explicit IndependentDraws(ModeDraws mode_draws)
        : mode_draws_{std::move(mode_draws)}
    {
    }

    auto append(const GeneticState& state) -> void
    {
        mode_draws_.for_each([&]<GeneticMode Mode>(auto& draws)
                             { draws.append(state.template get<Mode>()); });
    }

    template <GeneticMode Mode>
    [[nodiscard]] auto get() const noexcept -> const
        typename ModeDraws::template mode_value_type<Mode>&
    {
        return mode_draws_.template get<Mode>();
    }

   private:
    ModeDraws mode_draws_;
};

template <typename GeneticState, typename ModeDraws, typename JointT>
class JointDraws
{
   public:
    JointDraws(ModeDraws mode_draws, JointT joint)
        : mode_draws_{std::move(mode_draws)}, joint_{std::move(joint)}
    {
    }

    auto append(const GeneticState& state) -> void
    {
        mode_draws_.for_each(
            [&]<GeneticMode Mode>(auto& draws)
            { draws.append(state.mode_values().template get<Mode>()); });
        joint_.append(state.joint());
    }

    template <GeneticMode Mode>
    [[nodiscard]] auto get() const noexcept -> const
        typename ModeDraws::template mode_value_type<Mode>&
    {
        return mode_draws_.template get<Mode>();
    }

    [[nodiscard]] auto joint() const noexcept -> const JointT&
    {
        return joint_;
    }

   private:
    ModeDraws mode_draws_;
    JointT joint_;
};

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_DRAWS_H_
