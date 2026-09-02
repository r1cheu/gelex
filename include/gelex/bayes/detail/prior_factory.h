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

#ifndef GELEX_BAYES_DETAIL_PRIOR_FACTORY_H_
#define GELEX_BAYES_DETAIL_PRIOR_FACTORY_H_

#include <cmath>
#include <fmt/format.h>
#include <ranges>
#include <utility>
#include <vector>

#include "gelex/bayes/detail/calibration.h"
#include "gelex/bayes/detail/genetic_spec.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/variance_budget.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"
#include "gelex/infra/var.h"

namespace gelex
{

namespace detail
{

template <UpdatePolicy Policy, typename T, typename HyperPrior>
auto make_parameter(T initial, HyperPrior hyperprior)
{
    if constexpr (Policy == UpdatePolicy::Fixed)
    {
        return FixedParameter<T>{.initial = std::move(initial)};
    }
    else
    {
        return SampledParameter<T, HyperPrior>{
            .initial = std::move(initial), .hyperprior = std::move(hyperprior)};
    }
}

constexpr auto initial_activity(const ScaledMixture& spec) -> double
{
    auto activity = 0.0;
    for (const auto [probability, scale] :
         std::views::zip(spec.probabilities(), spec.scales()))
    {
        activity += probability * scale;
    }
    return activity;
}

template <GeneticMode Mode>
    requires(Mode == GeneticMode::A || Mode == GeneticMode::D)
constexpr auto initial_activity(const JointSpikeSlab& spec) -> double
{
    const auto& probabilities = spec.probabilities();
    if constexpr (Mode == GeneticMode::A)
    {
        return probabilities.at(1) + probabilities.at(3);
    }
    else
    {
        return probabilities.at(2) + probabilities.at(3);
    }
}

template <GeneticModeSet Modes, VarianceLayout Kind>
auto make_genetic_prior(
    GaussianFamily<Kind> /*family*/,
    const Gaussian& /*genetic_spec*/,
    const MarkerVarianceCalibrator& calibrator)
{
    return generate_mode_values<Modes>(
        [&]<GeneticMode Mode>() -> GaussianPrior<Kind>
        { return {.variance = calibrator.calibrate(Mode, 1.0)}; });
}

template <
    GeneticModeSet Modes,
    VarianceLayout Kind,
    UpdatePolicy ProbabilityUpdate>
auto make_genetic_prior(
    SpikeSlabFamily<Kind, ProbabilityUpdate> /*family*/,
    const genetic_spec_t<Modes, SpikeSlabFamily<Kind, ProbabilityUpdate>>&
        genetic_spec,
    const MarkerVarianceCalibrator& calibrator)
{
    return transform_mode_values(
        genetic_spec,
        [&]<GeneticMode Mode>(
            const SpikeSlab& spec) -> SpikeSlabPrior<Kind, ProbabilityUpdate>
        {
            return {
                .variance = calibrator.calibrate(Mode, spec.probability()),
                .probability = make_parameter<ProbabilityUpdate>(
                    spec.probability(), BetaHyperPrior{})};
        });
}

template <GeneticModeSet Modes, UpdatePolicy ProbabilitiesUpdate>
auto make_genetic_prior(
    ScaledMixtureFamily<ProbabilitiesUpdate> /*family*/,
    const genetic_spec_t<Modes, ScaledMixtureFamily<ProbabilitiesUpdate>>&
        genetic_spec,
    const MarkerVarianceCalibrator& calibrator)
{
    return transform_mode_values(
        genetic_spec,
        [&]<GeneticMode Mode>(const ScaledMixture& spec)
            -> ScaledMixturePrior<
                ScaledMixture::class_count,
                ProbabilitiesUpdate>
        {
            return {
                .variance = calibrator.calibrate(Mode, initial_activity(spec)),
                .probabilities = make_parameter<ProbabilitiesUpdate>(
                    spec.probabilities(),
                    DirichletHyperPrior<ScaledMixture::class_count>{}),
                .scales = spec.scales()};
        });
}

template <
    GeneticModeSet Modes,
    UpdatePolicy ProbabilitiesUpdate,
    HalfNormalAsymmetry Axis>
    requires(Modes == (GeneticMode::A | GeneticMode::D))
auto make_genetic_prior(
    JointSpikeSlabFamily<ProbabilitiesUpdate, Axis> /*family*/,
    const genetic_spec_t<
        Modes,
        JointSpikeSlabFamily<ProbabilitiesUpdate, Axis>>& genetic_spec,
    const MarkerVarianceCalibrator& calibrator)
{
    const auto& joint_spec = genetic_spec.joint();
    auto mode_priors = transform_mode_values(
        genetic_spec.mode_values(),
        [&]<GeneticMode Mode>(const auto& mode_spec)
        {
            auto variance = calibrator.calibrate(
                Mode, initial_activity<Mode>(joint_spec));
            if constexpr (Mode == GeneticMode::A)
            {
                return GaussianPrior<VarianceLayout::Pooled>{
                    .variance = std::move(variance)};
            }
            else if constexpr (Axis == HalfNormalAsymmetry::Count)
            {
                return HalfNormalPrior<Axis>{
                    .variance = std::move(variance),
                    .positive_probability
                    = make_parameter<UpdatePolicy::Sampled>(
                        mode_spec.positive_probability(), BetaHyperPrior{})};
            }
            else
            {
                return HalfNormalPrior<Axis>{.variance = std::move(variance)};
            }
        });

    using JointPrior
        = JointSpikeSlabPrior<JointSpikeSlab::class_count, ProbabilitiesUpdate>;
    return JointModeValues{
        std::move(mode_priors),
        JointPrior{
            .probabilities = make_parameter<ProbabilitiesUpdate>(
                joint_spec.probabilities(),
                DirichletHyperPrior<JointPrior::class_count>{})}};
}

inline auto random_projection_variance(const bayes::RandomDesign& design)
    -> double
{
    const double variance = matvar(design.X(), VarNormType::Population).sum();
    if (!std::isfinite(variance) || variance <= 0.0)
    {
        throw GelexException(
            fmt::format(
                "random design '{}' projection variance must be finite and "
                "positive, got {}",
                design.name(),
                variance));
    }
    return variance;
}

inline auto make_random_prior(
    const BayesModel& model,
    const VarianceBudget& budget,
    double phenotype_variance) -> std::vector<VarianceParameter>
{
    const auto designs = model.random();
    const double share = budget.random();
    if (designs.empty())
    {
        if (share != 0.0)
        {
            throw GelexException(
                "random variance share must be zero when the model has no "
                "random designs");
        }
        return {};
    }
    if (share <= 0.0)
    {
        throw GelexException(
            "random variance share must be positive when the model has "
            "random designs");
    }

    const double block_target
        = phenotype_variance * share / static_cast<double>(designs.size());
    std::vector<VarianceParameter> parameters;
    parameters.reserve(designs.size());
    for (const auto& design : designs)
    {
        const double initial
            = block_target / random_projection_variance(design);
        if (!std::isfinite(initial) || initial <= 0.0)
        {
            throw GelexException(
                fmt::format(
                    "random coefficient variance for design '{}' must be "
                    "finite and positive, got {}",
                    design.name(),
                    initial));
        }
        parameters.push_back(make_mean_calibrated_variance_parameter(initial));
    }
    return parameters;
}

}  // namespace detail

}  // namespace gelex

#endif  // GELEX_BAYES_DETAIL_PRIOR_FACTORY_H_
