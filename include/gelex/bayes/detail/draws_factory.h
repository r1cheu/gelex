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

#ifndef GELEX_BAYES_DETAIL_DRAWS_FACTORY_H_
#define GELEX_BAYES_DETAIL_DRAWS_FACTORY_H_

#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <string>
#include <utility>

#include "gelex/bayes/genetic/detail/draws_support.h"
#include "gelex/bayes/genetic/draws.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_writer.h"

namespace gelex::detail
{

template <VarianceLayout Kind>
[[nodiscard]] auto make_draws(
    const GaussianPrior<Kind>& /*prior*/,
    GeneticDrawsBuilder& builder) -> GaussianDraws<Kind>
{
    return {.variance = make_marker_variance_draw<Kind>(builder)};
}

[[nodiscard]] inline auto make_draws(
    const HalfNormalPrior& /*prior*/,
    GeneticDrawsBuilder& builder) -> HalfNormalDraws
{
    return {
        .variance = builder.scalar("variance"),
        .probit_coefficients = builder.vector("probit_coefficients", 2)};
}

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_draws(
    const SpikeSlabPrior<Kind, WeightUpdate>& /*prior*/,
    GeneticDrawsBuilder& builder) -> SpikeSlabDraws<Kind, WeightUpdate>
{
    return {
        .variance = make_marker_variance_draw<Kind>(builder),
        .assignment = builder.category<2>("assignment", builder.marker_count()),
        .probability
        = make_probability_draw<WeightUpdate>(builder, "probability")};
}

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_draws(
    const ScaledMixturePrior<ClassCount, WeightUpdate>& /*prior*/,
    GeneticDrawsBuilder& builder)
    -> ScaledMixtureDraws<ClassCount, WeightUpdate>
{
    return {
        .variance = builder.scalar("variance"),
        .assignment
        = builder.category<ClassCount>("assignment", builder.marker_count()),
        .probabilities
        = make_probabilities_draw<WeightUpdate, ClassCount>(builder),
        .component_explained_variance = make_component_explained_variance_draw<
            ScaledMixtureState<ClassCount>::component_count>(builder)};
}

// Rows follow the fitted column layout of JointSpikeSlabState: A in A-only,
// A in AD, D in D-only, D in AD.
template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_draws(
    const JointSpikeSlabPrior<ClassCount, WeightUpdate>& /*prior*/,
    GeneticDrawsBuilder& builder)
    -> JointSpikeSlabDraws<ClassCount, WeightUpdate>
{
    return {
        .assignment
        = builder.category<ClassCount>("assignment", builder.marker_count()),
        .probabilities
        = make_probabilities_draw<WeightUpdate, ClassCount>(builder),
        .component_explained_variance = make_component_explained_variance_draw<
            JointSpikeSlabState<ClassCount>::component_count>(builder)};
}

template <typename Prior>
[[nodiscard]] auto make_family_draws(
    const Prior& prior,
    BinaryWriter& writer,
    std::string prefix,
    GeneticDrawsDimensions dimensions)
{
    GeneticDrawsBuilder builder{writer, std::move(prefix), dimensions};
    return make_draws(prior, builder);
}

template <GeneticModeSet Modes>
[[nodiscard]] auto make_coefficient_draws(
    BinaryWriter& writer,
    const GeneticDrawsDimensions& dimensions)
{
    auto draws = generate_mode_values<Modes>(
        [&]<GeneticMode Mode>()
        {
            GeneticDrawsBuilder builder{
                writer, fmt::format("genetic/{}", Mode), dimensions};
            return builder.vector("coefficients", dimensions.marker_count);
        });
    return GeneticCoefficientDraws<Modes>{std::move(draws)};
}

template <GeneticModeSet Modes, typename... Priors>
[[nodiscard]] auto make_mode_family_draws(
    const ModeValues<Modes, Priors...>& prior,
    BinaryWriter& writer,
    GeneticDrawsDimensions dimensions)
{
    return transform_mode_values(
        prior,
        [&]<GeneticMode Mode>(const auto& mode_prior)
        {
            return make_family_draws(
                mode_prior,
                writer,
                fmt::format("genetic/{}", Mode),
                dimensions);
        });
}

template <GeneticModeSet Modes, typename... Priors>
[[nodiscard]] auto make_draws(
    const ModeValues<Modes, Priors...>& prior,
    const bayes::GeneticDesign& design,
    BinaryWriter& writer,
    std::uint64_t draw_count)
{
    const auto dimensions = GeneticDrawsDimensions{
        .marker_count = design.cols(), .draw_count = draw_count};
    auto coefficient_draws = make_coefficient_draws<Modes>(writer, dimensions);
    auto family_draws = make_mode_family_draws(prior, writer, dimensions);
    return IndependentGeneticDraws<
        decltype(coefficient_draws),
        decltype(family_draws)>{
        std::move(coefficient_draws), std::move(family_draws)};
}

template <typename ModeValuesType, typename JointPrior>
[[nodiscard]] auto make_draws(
    const JointModeValues<ModeValuesType, JointPrior>& prior,
    const bayes::GeneticDesign& design,
    BinaryWriter& writer,
    std::uint64_t draw_count)
{
    const auto dimensions = GeneticDrawsDimensions{
        .marker_count = design.cols(), .draw_count = draw_count};
    auto coefficient_draws
        = make_coefficient_draws<ModeValuesType::modes>(writer, dimensions);
    auto family_draws
        = make_mode_family_draws(prior.mode_values(), writer, dimensions);
    GeneticDrawsBuilder joint_builder{writer, "genetic/joint", dimensions};
    auto joint_draws = make_draws(prior.joint(), joint_builder);
    return JointGeneticDraws<
        decltype(coefficient_draws),
        decltype(family_draws),
        decltype(joint_draws)>{
        std::move(coefficient_draws),
        std::move(family_draws),
        std::move(joint_draws)};
}

template <typename Prior>
using genetic_draws_t = decltype(make_draws(
    std::declval<const Prior&>(),
    std::declval<const bayes::GeneticDesign&>(),
    std::declval<BinaryWriter&>(),
    std::declval<std::uint64_t>()));

}  // namespace gelex::detail

#endif  // GELEX_BAYES_DETAIL_DRAWS_FACTORY_H_
