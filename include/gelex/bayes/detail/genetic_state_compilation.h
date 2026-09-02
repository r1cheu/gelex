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
            prior.variance, dimensions.marker_count),
        .fitted_values = Eigen::VectorXd::Zero(dimensions.individual_count)};
}

inline auto make_parameter_state(
    const HalfNormalPrior<HalfNormalAsymmetry::Count>& prior,
    GeneticStateDimensions dimensions)
    -> HalfNormalState<HalfNormalAsymmetry::Count>
{
    return {
        .variance = prior.variance.initial_value(),
        .assignment
        = Eigen::VectorX<std::uint8_t>::Zero(dimensions.marker_count),
        .positive_probability = prior.positive_probability.initial,
        .fitted_values = Eigen::VectorXd::Zero(dimensions.individual_count)};
}

inline auto make_parameter_state(
    const HalfNormalPrior<HalfNormalAsymmetry::Magnitude>& prior,
    GeneticStateDimensions dimensions)
    -> HalfNormalState<HalfNormalAsymmetry::Magnitude>
{
    const double initial_variance = prior.variance.initial_value();
    return {
        .variances = {initial_variance, initial_variance},
        .assignment
        = Eigen::VectorX<std::uint8_t>::Zero(dimensions.marker_count),
        .fitted_values = Eigen::VectorXd::Zero(dimensions.individual_count)};
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
        .probability = prior.probability.initial,
        .fitted_values = Eigen::VectorXd::Zero(dimensions.individual_count)};
}

template <std::size_t ClassCount, UpdatePolicy ProbabilitiesUpdate>
auto make_parameter_state(
    const ScaledMixturePrior<ClassCount, ProbabilitiesUpdate>& prior,
    GeneticStateDimensions dimensions) -> ScaledMixtureState<ClassCount>
{
    auto state = ScaledMixtureState<ClassCount>{
        .variance = prior.variance.initial_value(),
        .assignment
        = Eigen::VectorX<std::uint8_t>::Zero(dimensions.marker_count),
        .probabilities = prior.probabilities.initial,
    };
    state.fitted_values.setZero(
        dimensions.individual_count,
        static_cast<Eigen::Index>(state.component_count));
    return state;
}

template <std::size_t ClassCount, UpdatePolicy ProbabilitiesUpdate>
auto make_parameter_state(
    const JointSpikeSlabPrior<ClassCount, ProbabilitiesUpdate>& prior,
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

template <typename Prior>
auto make_mode_state(const Prior& prior, GeneticStateDimensions dimensions)
{
    auto family_state = make_parameter_state(prior, dimensions);
    using State = GeneticModeState<decltype(family_state)>;
    return State{
        .coefficients = Eigen::VectorXd::Zero(dimensions.marker_count),
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

template <GeneticModeSet Modes, typename... Priors>
auto make_genetic_state(
    const IndependentTopology<Modes, Priors...>& prior,
    const bayes::GeneticDesign& design)
{
    validate_state_design<Modes>(design);
    const auto dimensions = GeneticStateDimensions{
        .marker_count = design.cols(), .individual_count = design.rows()};
    return transform_mode_values(
        prior,
        [&]<GeneticMode /*Mode*/>(const auto& mode_prior)
        { return make_mode_state(mode_prior, dimensions); });
}

template <typename ModeTopology, typename JointPrior>
auto make_genetic_state(
    const JointTopology<ModeTopology, JointPrior>& prior,
    const bayes::GeneticDesign& design)
{
    validate_state_design<JointTopology<ModeTopology, JointPrior>::modes>(
        design);
    const auto dimensions = GeneticStateDimensions{
        .marker_count = design.cols(), .individual_count = design.rows()};
    auto mode_states = transform_mode_values(
        prior.mode_values(),
        [&]<GeneticMode /*Mode*/>(const auto& mode_prior)
        { return make_mode_state(mode_prior, dimensions); });
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
