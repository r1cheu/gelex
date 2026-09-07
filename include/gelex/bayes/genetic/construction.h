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

#ifndef GELEX_BAYES_GENETIC_CONSTRUCTION_H_
#define GELEX_BAYES_GENETIC_CONSTRUCTION_H_

// Composes per-method prior, state, and draw factories over genetic modes.

#include <Eigen/Core>
#include <cstddef>
#include <cstdint>
#include <utility>

#include "gelex/bayes/genetic/draws.h"
#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/joint_spike_slab.h"
#include "gelex/bayes/genetic/scaled_mixture.h"
#include "gelex/bayes/genetic/spike_slab.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/variance/detail/calibration.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_writer.h"

namespace gelex::detail
{

// ---- prior

template <GeneticModeSet Modes, typename... Specs>
auto make_prior(
    const ModeValues<Modes, Specs...>& specs,
    const MarkerVarianceCalibrator& calibrator)
{
    return transform_mode_values(
        specs,
        [&]<GeneticMode Mode>(const auto& spec)
        { return make_prior<Mode>(spec, calibrator); });
}

// ---- state

template <GeneticModeSet Modes>
auto validate_genetic_design(const bayes::GeneticDesign& design) -> void
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
    validate_genetic_design<Modes>(design);
    const auto dimensions = GeneticDimensions{
        .individual = static_cast<std::size_t>(design.rows()),
        .marker = static_cast<std::size_t>(design.cols())};
    return transform_mode_values(
        prior,
        [&]<GeneticMode /*Mode*/>(const auto& mode_prior)
        { return make_state(mode_prior, dimensions); });
}

template <typename ModeValuesType, typename JointPrior>
auto make_state(
    const JointModeValues<ModeValuesType, JointPrior>& prior,
    const bayes::GeneticDesign& design)
{
    auto mode_states = make_state(prior.mode_values(), design);
    const auto dimensions = GeneticDimensions{
        .individual = static_cast<std::size_t>(design.rows()),
        .marker = static_cast<std::size_t>(design.cols())};
    return JointModeValues{
        std::move(mode_states), make_state(prior.joint(), dimensions)};
}

template <typename Prior>
using genetic_state_t = decltype(make_state(
    std::declval<const Prior&>(),
    std::declval<const bayes::GeneticDesign&>()));

// ---- draws

template <GeneticModeSet Modes, typename... Priors>
[[nodiscard]] auto make_draws(
    const ModeValues<Modes, Priors...>& prior,
    const bayes::GeneticDesign& design,
    BinaryWriter& writer,
    std::uint64_t draw_count)
{
    validate_genetic_design<Modes>(design);
    const GeneticDimensions dimensions{
        .individual = static_cast<std::size_t>(design.rows()),
        .marker = static_cast<std::size_t>(design.cols())};
    return transform_mode_values(
        prior,
        [&]<GeneticMode Mode>(const auto& mode_prior)
        {
            return make_draws(
                mode_prior, writer, genetic_id<Mode>, draw_count, dimensions);
        });
}

template <typename ModeValuesType, typename JointPrior>
[[nodiscard]] auto make_draws(
    const JointModeValues<ModeValuesType, JointPrior>& prior,
    const bayes::GeneticDesign& design,
    BinaryWriter& writer,
    std::uint64_t draw_count)
{
    auto mode_draws
        = make_draws(prior.mode_values(), design, writer, draw_count);
    const GeneticDimensions dimensions{
        .individual = static_cast<std::size_t>(design.rows()),
        .marker = static_cast<std::size_t>(design.cols())};
    auto joint_draws = make_draws(
        prior.joint(), writer, joint_genetic_id, draw_count, dimensions);
    return JointModeValues{std::move(mode_draws), std::move(joint_draws)};
}

template <typename Prior>
using genetic_draws_t = decltype(make_draws(
    std::declval<const Prior&>(),
    std::declval<const bayes::GeneticDesign&>(),
    std::declval<BinaryWriter&>(),
    std::declval<std::uint64_t>()));

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_CONSTRUCTION_H_
