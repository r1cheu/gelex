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

#ifndef GELEX_BAYES_PRIOR_COMPILATION_H_
#define GELEX_BAYES_PRIOR_COMPILATION_H_

#include <cstddef>
#include <optional>
#include <ranges>
#include <utility>

#include "gelex/bayes/detail/calibration.h"
#include "gelex/bayes/genetic/independent_topology.h"
#include "gelex/bayes/genetic/joint_topology.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/semantic_method.h"
#include "gelex/bayes/spec.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{

namespace detail
{

template <typename Value, typename Hyper>
auto compile_updatable(const MaybeSampled<Value>& spec, Hyper hyper)
    -> Updatable<Value, Hyper>
{
    return {
        .initial = spec.initial,
        .prior = spec.update == UpdatePolicy::Sampled
                     ? std::optional<Hyper>{hyper}
                     : std::nullopt};
}

constexpr auto initial_activity(const ScaledMixture& spec) -> double
{
    auto activity = 0.0;
    for (const auto [probability, scale] :
         std::views::zip(spec.probabilities.initial, spec.scales))
    {
        activity += probability * scale;
    }
    return activity;
}

template <GeneticModeSet Modes, Variance Kind>
auto compile_genetic_prior(
    GaussianMethod<Kind> /*method*/,
    const NoParameters& /*parameters*/,
    const MarkerVarianceCalibrator& calibrator)
    -> IndependentTopology<Modes, GaussianPrior<Kind>>
{
    return generate_mode_values<Modes>(
        [&](GeneticMode mode) -> GaussianPrior<Kind>
        { return {.variance = calibrator.calibrate(mode, 1.0)}; });
}

template <GeneticModeSet Modes, Variance Kind>
auto compile_genetic_prior(
    SpikeSlabMethod<Kind> /*method*/,
    const IndependentTopology<Modes, SpikeSlab>& parameters,
    const MarkerVarianceCalibrator& calibrator)
    -> IndependentTopology<Modes, SpikeSlabPrior<Kind>>
{
    return transform_mode_values(
        parameters,
        [&](GeneticMode mode, const SpikeSlab& spec) -> SpikeSlabPrior<Kind>
        {
            return {
                .variance
                = calibrator.calibrate(mode, spec.probability.initial),
                .probability
                = compile_updatable(spec.probability, BetaPrior{})};
        });
}

template <GeneticModeSet Modes>
auto compile_genetic_prior(
    ScaledMixtureMethod /*method*/,
    const IndependentTopology<Modes, ScaledMixture>& parameters,
    const MarkerVarianceCalibrator& calibrator)
    -> IndependentTopology<Modes, ScaledMixturePrior>
{
    return transform_mode_values(
        parameters,
        [&](GeneticMode mode, const ScaledMixture& spec) -> ScaledMixturePrior
        {
            return {
                .variance = calibrator.calibrate(mode, initial_activity(spec)),
                .probabilities
                = compile_updatable(spec.probabilities, DirichletPrior<5>{}),
                .scales = spec.scales};
        });
}

template <GeneticModeSet Modes>
    requires(Modes == (GeneticMode::A | GeneticMode::D))
auto compile_genetic_prior(
    JointSpikeSlabMethod /*method*/,
    const JointSpikeSlab& parameters,
    const MarkerVarianceCalibrator& calibrator)
    -> JointTopology<GaussianPrior<Variance::Pooled>, JointSpikeSlabPrior>
{
    constexpr std::size_t A_ONLY = 1;
    constexpr std::size_t D_ONLY = 2;
    constexpr std::size_t BOTH = 3;

    const auto& probabilities = parameters.probabilities.initial;
    auto mode_values = generate_mode_values<Modes>(
        [&](GeneticMode mode) -> GaussianPrior<Variance::Pooled>
        {
            const double activity
                = mode == GeneticMode::A
                      ? probabilities[A_ONLY] + probabilities[BOTH]
                      : probabilities[D_ONLY] + probabilities[BOTH];
            return {.variance = calibrator.calibrate(mode, activity)};
        });

    return JointTopology{
        std::move(mode_values),
        JointSpikeSlabPrior{
            .probabilities
            = compile_updatable(parameters.probabilities, DirichletPrior<4>{}),
            .positive_probability
            = compile_updatable(parameters.positive_probability, BetaPrior{})}};
}

}  // namespace detail

template <typename SemanticMethod, GeneticModeSet Modes>
auto compile_genetic_prior(
    const BayesRecipe<SemanticMethod, Modes>& recipe,
    const BayesModel& model)
{
    const detail::MarkerVarianceCalibrator calibrator{model, recipe.variance()};
    return detail::compile_genetic_prior<Modes>(
        SemanticMethod{}, recipe.parameters(), calibrator);
}

template <typename Recipe>
using genetic_prior_t = decltype(compile_genetic_prior(
    std::declval<const Recipe&>(),
    std::declval<const BayesModel&>()));

}  // namespace gelex

#endif  // GELEX_BAYES_PRIOR_COMPILATION_H_
