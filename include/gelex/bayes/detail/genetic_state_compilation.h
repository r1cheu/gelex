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

#ifndef GELEX_BAYES_DETAIL_GENETIC_STATE_COMPILATION_H_
#define GELEX_BAYES_DETAIL_GENETIC_STATE_COMPILATION_H_

#include <Eigen/Core>
#include <cstdint>
#include <utility>

#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/independent_topology.h"
#include "gelex/bayes/genetic/joint_topology.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/exception.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::detail
{

struct GeneticStateDimensions
{
    Eigen::Index marker_count;
    Eigen::Index individual_count;
};

template <VarianceLayout Kind>
auto initial_marker_variance(
    const VarianceParameter& parameter,
    Eigen::Index marker_count) -> marker_variance_state_t<Kind>
{
    if constexpr (Kind == VarianceLayout::Pooled)
    {
        return parameter.initial_value();
    }
    else
    {
        return Eigen::VectorXd::Constant(
            marker_count, parameter.initial_value());
    }
}

template <VarianceLayout Kind>
auto make_parameter_state(
    const GaussianPrior<Kind>& prior,
    GeneticStateDimensions dimensions) -> GaussianState<Kind>
{
    return {
        .variance = initial_marker_variance<Kind>(
            prior.variance, dimensions.marker_count)};
}

template <VarianceLayout Kind, UpdatePolicy ProbabilityUpdate>
auto make_parameter_state(
    const SpikeSlabPrior<Kind, ProbabilityUpdate>& prior,
    GeneticStateDimensions dimensions) -> SpikeSlabState<Kind>
{
    return {
        .variance = initial_marker_variance<Kind>(
            prior.variance, dimensions.marker_count),
        .assignment
        = Eigen::VectorX<std::uint8_t>::Zero(dimensions.marker_count),
        .probability = prior.probability.initial};
}

template <UpdatePolicy ProbabilitiesUpdate>
auto make_parameter_state(
    const ScaledMixturePrior<ProbabilitiesUpdate>& prior,
    GeneticStateDimensions dimensions) -> ScaledMixtureState
{
    return {
        .variance = prior.variance.initial_value(),
        .assignment
        = Eigen::VectorX<std::uint8_t>::Zero(dimensions.marker_count),
        .probabilities = prior.probabilities.initial,
    };
}

template <
    UpdatePolicy ProbabilitiesUpdate,
    UpdatePolicy PositiveProbabilityUpdate>
auto make_parameter_state(
    const JointSpikeSlabPrior<ProbabilitiesUpdate, PositiveProbabilityUpdate>&
        prior,
    GeneticStateDimensions dimensions) -> JointSpikeSlabState
{
    return {
        .assignment
        = Eigen::VectorX<std::uint8_t>::Zero(dimensions.marker_count),
        .probabilities = prior.probabilities.initial,

        .dominance_sign
        = Eigen::VectorX<std::uint8_t>::Zero(dimensions.marker_count),
        .positive_probability = prior.positive_probability.initial};
}

template <ComponentLayout Layout, typename Prior>
auto make_mode_state(const Prior& prior, GeneticStateDimensions dimensions)
{
    auto family_state = make_parameter_state(prior, dimensions);
    using State = GeneticModeState<decltype(family_state), Layout>;
    using FittedMatrix = typename State::ComponentFittedMatrix;
    return State{
        .coefficients = Eigen::VectorXd::Zero(dimensions.marker_count),
        .component_fitted_values = FittedMatrix::Zero(
            dimensions.individual_count,
            static_cast<Eigen::Index>(Layout::component_count)),
        .family_state = std::move(family_state)};
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

template <GeneticModeSet Modes, typename Prior>
auto make_genetic_state(
    const IndependentTopology<Modes, Prior>& prior,
    const bayes::GeneticDesign& design)
{
    using Layout = typename Prior::component_layout;
    validate_state_design<Modes>(design);
    const auto dimensions = GeneticStateDimensions{
        .marker_count = design.cols(), .individual_count = design.rows()};
    return transform_mode_values(
        prior,
        [&](GeneticMode /*mode*/, const Prior& mode_prior)
        { return make_mode_state<Layout>(mode_prior, dimensions); });
}

template <typename ModePrior, typename JointPrior>
auto make_genetic_state(
    const JointTopology<ModePrior, JointPrior>& prior,
    const bayes::GeneticDesign& design)
{
    using Layout = typename JointPrior::component_layout;
    validate_state_design<JointTopology<ModePrior, JointPrior>::modes>(design);
    const auto dimensions = GeneticStateDimensions{
        .marker_count = design.cols(), .individual_count = design.rows()};
    auto mode_states = transform_mode_values(
        prior.mode_values(),
        [&](GeneticMode /*mode*/, const ModePrior& mode_prior)
        { return make_mode_state<Layout>(mode_prior, dimensions); });
    return JointTopology{
        std::move(mode_states),
        make_parameter_state(prior.joint(), dimensions)};
}

template <typename Prior>
using genetic_state_t = decltype(make_genetic_state(
    std::declval<const Prior&>(),
    std::declval<const bayes::GeneticDesign&>()));

}  // namespace gelex::detail

#endif  // GELEX_BAYES_DETAIL_GENETIC_STATE_COMPILATION_H_
