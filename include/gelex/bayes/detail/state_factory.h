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

#ifndef GELEX_BAYES_DETAIL_STATE_FACTORY_H_
#define GELEX_BAYES_DETAIL_STATE_FACTORY_H_

#include <Eigen/Core>
#include <cstdint>
#include <utility>

#include "gelex/bayes/genetic/detail/state_support.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

namespace gelex::detail
{

template <VarianceLayout Kind>
auto make_state(
    const GaussianPrior<Kind>& prior,
    GeneticStateDimensions dimensions) -> GaussianState<Kind>
{
    return {
        .variance = initial_marker_variance<Kind>(
            prior.variance, dimensions.marker_count),
        .fitted_values = Eigen::VectorXd::Zero(dimensions.individual_count)};
}

inline auto make_state(
    const HalfNormalPrior& prior,
    GeneticStateDimensions dimensions) -> HalfNormalState
{
    return {
        .variance = prior.variance.initial,
        .probit_coefficients = Eigen::Vector2d::Zero(),
        .fitted_values = Eigen::VectorXd::Zero(dimensions.individual_count)};
}

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
auto make_state(
    const SpikeSlabPrior<Kind, WeightUpdate>& prior,
    GeneticStateDimensions dimensions) -> SpikeSlabState<Kind>
{
    return {
        .variance = initial_marker_variance<Kind>(
            prior.variance, dimensions.marker_count),
        .assignment
        = Eigen::VectorX<std::uint8_t>::Zero(dimensions.marker_count),
        .probability = prior.probability.initial,
        .fitted_values = Eigen::VectorXd::Zero(dimensions.individual_count)};
}

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
auto make_state(
    const ScaledMixturePrior<ClassCount, WeightUpdate>& prior,
    GeneticStateDimensions dimensions) -> ScaledMixtureState<ClassCount>
{
    auto state = ScaledMixtureState<ClassCount>{
        .variance = prior.variance.initial,
        .assignment
        = Eigen::VectorX<std::uint8_t>::Zero(dimensions.marker_count),
        .probabilities = prior.probabilities.initial,
    };
    state.fitted_values.setZero(
        dimensions.individual_count,
        static_cast<Eigen::Index>(state.component_count));
    return state;
}

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
auto make_state(
    const JointSpikeSlabPrior<ClassCount, WeightUpdate>& prior,
    GeneticStateDimensions dimensions) -> JointSpikeSlabState<ClassCount>
{
    auto state = JointSpikeSlabState<ClassCount>{
        .assignment
        = Eigen::VectorX<std::uint8_t>::Zero(dimensions.marker_count),
        .probabilities = prior.probabilities.initial,
        .fitted_values = {},
    };
    state.fitted_values.setZero(
        dimensions.individual_count,
        static_cast<Eigen::Index>(state.component_count));
    return state;
}

template <GeneticModeSet Modes>
auto validate_state_design(const bayes::GeneticDesign& design) -> void
{
    for (const auto mode : Modes.each())
    {
        if (!design.contains(mode))
        {
            throw GelexException(
                "genetic design does not contain every mode required by the "
                "prior");
        }
    }
}

template <GeneticModeSet Modes, typename... Priors>
auto make_state(
    const ModeValues<Modes, Priors...>& prior,
    const bayes::GeneticDesign& design)
{
    validate_state_design<Modes>(design);
    const auto dimensions = GeneticStateDimensions{
        .marker_count = design.cols(), .individual_count = design.rows()};
    return transform_mode_values(
        prior,
        [&]<GeneticMode /*Mode*/>(const auto& mode_prior)
        {
            auto family_state = make_state(mode_prior, dimensions);
            return GeneticModeState<decltype(family_state)>{
                .coefficients = Eigen::VectorXd::Zero(dimensions.marker_count),
                .family_state = std::move(family_state)};
        });
}

template <typename ModeValuesType, typename JointPrior>
auto make_state(
    const JointModeValues<ModeValuesType, JointPrior>& prior,
    const bayes::GeneticDesign& design)
{
    auto mode_states = make_state(prior.mode_values(), design);
    const auto dimensions = GeneticStateDimensions{
        .marker_count = design.cols(), .individual_count = design.rows()};
    return JointModeValues{
        std::move(mode_states), make_state(prior.joint(), dimensions)};
}

template <typename Prior>
using genetic_state_t = decltype(make_state(
    std::declval<const Prior&>(),
    std::declval<const bayes::GeneticDesign&>()));

}  // namespace gelex::detail

#endif  // GELEX_BAYES_DETAIL_STATE_FACTORY_H_
