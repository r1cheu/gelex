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

#include <Eigen/Core>
#include <cmath>
#include <fmt/format.h>
#include <ranges>
#include <utility>
#include <vector>

#include "gelex/bayes/detail/calibration.h"
#include "gelex/bayes/genetic/independent_topology.h"
#include "gelex/bayes/genetic/joint_topology.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/semantic_method.h"
#include "gelex/bayes/spec.h"
#include "gelex/exception.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{

namespace detail
{

template <UpdatePolicy Policy, typename T, typename HyperPrior>
auto compile_parameter(T initial, HyperPrior hyperprior)
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
auto compile_genetic_prior(
    GaussianMethod<Kind> /*method*/,
    const Gaussian& /*parameters*/,
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
auto compile_genetic_prior(
    SpikeSlabMethod<Kind, ProbabilityUpdate> /*method*/,
    const method_parameters_t<SpikeSlabMethod<Kind, ProbabilityUpdate>, Modes>&
        parameters,
    const MarkerVarianceCalibrator& calibrator)
{
    return transform_mode_values(
        parameters,
        [&]<GeneticMode Mode>(
            const SpikeSlab& spec) -> SpikeSlabPrior<Kind, ProbabilityUpdate>
        {
            return {
                .variance = calibrator.calibrate(Mode, spec.probability()),
                .probability = compile_parameter<ProbabilityUpdate>(
                    spec.probability(), BetaHyperPrior{})};
        });
}

template <GeneticModeSet Modes, UpdatePolicy ProbabilitiesUpdate>
auto compile_genetic_prior(
    ScaledMixtureMethod<ProbabilitiesUpdate> /*method*/,
    const method_parameters_t<ScaledMixtureMethod<ProbabilitiesUpdate>, Modes>&
        parameters,
    const MarkerVarianceCalibrator& calibrator)
{
    return transform_mode_values(
        parameters,
        [&]<GeneticMode Mode>(const ScaledMixture& spec)
            -> ScaledMixturePrior<
                ScaledMixture::class_count,
                ProbabilitiesUpdate>
        {
            return {
                .variance = calibrator.calibrate(Mode, initial_activity(spec)),
                .probabilities = compile_parameter<ProbabilitiesUpdate>(
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
auto compile_genetic_prior(
    JointSpikeSlabMethod<ProbabilitiesUpdate, Axis> /*method*/,
    const method_parameters_t<
        JointSpikeSlabMethod<ProbabilitiesUpdate, Axis>,
        Modes>& parameters,
    const MarkerVarianceCalibrator& calibrator)
{
    const auto& joint_parameters = parameters.joint();
    auto mode_priors = transform_mode_values(
        parameters.mode_values(),
        [&]<GeneticMode Mode>(const auto& mode_parameters)
        {
            auto variance = calibrator.calibrate(
                Mode, initial_activity<Mode>(joint_parameters));
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
                    = compile_parameter<UpdatePolicy::Sampled>(
                        mode_parameters.positive_probability(),
                        BetaHyperPrior{})};
            }
            else
            {
                return HalfNormalPrior<Axis>{.variance = std::move(variance)};
            }
        });

    using JointPrior
        = JointSpikeSlabPrior<JointSpikeSlab::class_count, ProbabilitiesUpdate>;
    return JointTopology{
        std::move(mode_priors),
        JointPrior{
            .probabilities = compile_parameter<ProbabilitiesUpdate>(
                joint_parameters.probabilities(),
                DirichletHyperPrior<JointPrior::class_count>{})}};
}

inline auto validate_phenotype_variance(const BayesModel& model) -> double
{
    const double variance = model.phenotype_variance();
    if (!std::isfinite(variance) || variance <= 0.0)
    {
        throw GelexException(
            fmt::format(
                "phenotype variance must be finite and positive, got {}",
                variance));
    }
    return variance;
}

inline auto random_projection_variance(const bayes::RandomDesign& design)
    -> double
{
    if (design.X.rows() == 0)
    {
        throw GelexException(
            fmt::format("random design '{}' has no rows", design.name));
    }

    double variance = 0.0;
    for (Eigen::Index column = 0; column < design.X.cols(); ++column)
    {
        const auto values = design.X.col(column);
        const double mean = values.mean();
        variance += (values.array() - mean).square().mean();
    }
    if (!std::isfinite(variance) || variance <= 0.0)
    {
        throw GelexException(
            fmt::format(
                "random design '{}' projection variance must be finite and "
                "positive, got {}",
                design.name,
                variance));
    }
    return variance;
}

inline auto compile_random_prior(
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
                    design.name,
                    initial));
        }
        parameters.push_back(make_mean_calibrated_variance_parameter(initial));
    }
    return parameters;
}

inline auto compile_residual_prior(
    const VarianceBudget& budget,
    double phenotype_variance) -> VarianceParameter
{
    const double initial = phenotype_variance * budget.residual();
    if (!std::isfinite(initial) || initial <= 0.0)
    {
        throw GelexException(
            fmt::format(
                "residual variance must be finite and positive, got {}",
                initial));
    }
    return make_mean_calibrated_variance_parameter(initial);
}

}  // namespace detail

template <typename SemanticMethod, GeneticModeSet Modes>
auto compile(
    const BayesRecipe<SemanticMethod, Modes>& recipe,
    const BayesModel& model)
{
    const double phenotype_variance
        = detail::validate_phenotype_variance(model);
    const detail::MarkerVarianceCalibrator calibrator{model, recipe.variance()};
    auto genetic = detail::compile_genetic_prior<Modes>(
        SemanticMethod{}, recipe.parameters(), calibrator);
    auto random = detail::compile_random_prior(
        model, recipe.variance(), phenotype_variance);
    auto residual
        = detail::compile_residual_prior(recipe.variance(), phenotype_variance);
    return detail::make_bayes_prior(
        std::move(random), std::move(genetic), residual);
}

template <typename Recipe>
using bayes_prior_t = decltype(compile(
    std::declval<const Recipe&>(),
    std::declval<const BayesModel&>()));

}  // namespace gelex

#endif  // GELEX_BAYES_PRIOR_COMPILATION_H_
