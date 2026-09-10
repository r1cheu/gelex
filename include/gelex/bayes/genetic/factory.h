// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_FACTORY_H_
#define GELEX_BAYES_GENETIC_FACTORY_H_

// Composes per-method prior, state, and draw factories over genetic modes.

#include <Eigen/Core>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "gelex/bayes/genetic/diagnostics_traits.h"
#include "gelex/bayes/genetic/draw_traits.h"
#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/joint_spike_slab.h"
#include "gelex/bayes/genetic/marker_effect_traits.h"
#include "gelex/bayes/genetic/scaled_mixture.h"
#include "gelex/bayes/genetic/spike_slab.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/serialization_ids.h"
#include "gelex/bayes/variance/calibration.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

namespace gelex
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

namespace detail
{
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

}  // namespace detail

template <GeneticModeSet Modes, typename... Priors>
auto make_state(
    const ModeValues<Modes, Priors...>& prior,
    const bayes::GeneticDesign& design)
{
    detail::validate_genetic_design<Modes>(design);
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

template <GeneticModeSet Modes, typename... States>
[[nodiscard]] auto make_draws(
    const ModeValues<Modes, States...>& state,
    DrawWriters writers,
    std::uint64_t draw_count)
{
    return transform_mode_values(
        state,
        [&]<GeneticMode Mode>(const auto& mode_state)
        {
            return make_draws(
                mode_state, writers, genetic_id<Mode>, draw_count);
        });
}

// The joint assignment zeroes inactive mode coefficients, so they are sparse.
template <typename ModeValuesType, typename JointState>
[[nodiscard]] auto make_draws(
    const JointModeValues<ModeValuesType, JointState>& state,
    DrawWriters writers,
    std::uint64_t draw_count)
{
    auto mode_draws = transform_mode_values(
        state.mode_values(),
        [&]<GeneticMode Mode>(const auto& mode_state)
        {
            return make_draws<CoefficientLayout::Sparse>(
                mode_state, writers, genetic_id<Mode>, draw_count);
        });
    auto joint_draws
        = make_draws(state.joint(), writers, joint_genetic_id, draw_count);
    return JointModeValues{std::move(mode_draws), std::move(joint_draws)};
}

template <typename State>
using genetic_draws_t = decltype(make_draws(
    std::declval<const State&>(),
    std::declval<DrawWriters>(),
    std::declval<std::uint64_t>()));

// ---- diagnostics

template <GeneticModeSet Modes, typename... Draws>
[[nodiscard]] auto make_diagnostics(
    std::type_identity<ModeValues<Modes, Draws...>> /*draws*/,
    DrawReaders readers,
    double prob)
{
    return generate_mode_values<Modes>(
        [&]<GeneticMode Mode>()
        {
            using draws_type =
                typename ModeValues<Modes, Draws...>::template mode_value_type<
                    Mode>;
            return make_diagnostics(
                std::type_identity<draws_type>{},
                readers,
                genetic_id<Mode>,
                prob);
        });
}

template <typename ModeDraws, typename JointDraws>
[[nodiscard]] auto make_diagnostics(
    std::type_identity<JointModeValues<ModeDraws, JointDraws>> /*draws*/,
    DrawReaders readers,
    double prob)
{
    auto mode_diagnostics
        = make_diagnostics(std::type_identity<ModeDraws>{}, readers, prob);
    auto joint_diagnostics = make_diagnostics(
        std::type_identity<JointDraws>{}, readers, joint_genetic_id, prob);
    return JointModeValues{
        std::move(mode_diagnostics), std::move(joint_diagnostics)};
}

template <GeneticModeSet Modes, typename... Diagnostics>
auto append_diagnostic_entries(
    std::vector<DiagnosticEntry>& out,
    const ModeValues<Modes, Diagnostics...>& diagnostics) -> void
{
    diagnostics.for_each(
        [&]<GeneticMode Mode>(const auto& mode_diagnostics)
        {
            append_diagnostic_entries(out, mode_diagnostics, genetic_id<Mode>);
        });
}

template <typename ModeDiagnostics, typename JointDiagnostics>
auto append_diagnostic_entries(
    std::vector<DiagnosticEntry>& out,
    const JointModeValues<ModeDiagnostics, JointDiagnostics>& diagnostics)
    -> void
{
    append_diagnostic_entries(out, diagnostics.mode_values());
    append_diagnostic_entries(out, diagnostics.joint(), joint_genetic_id);
}

// ---- marker effects

template <GeneticModeSet Modes, typename... Draws>
[[nodiscard]] auto make_marker_effects(
    std::type_identity<ModeValues<Modes, Draws...>> /*draws*/,
    DrawReaders readers)
{
    return generate_mode_values<Modes>(
        [&]<GeneticMode Mode>()
        {
            using draws_type =
                typename ModeValues<Modes, Draws...>::template mode_value_type<
                    Mode>;
            return make_marker_effects(
                std::type_identity<draws_type>{}, readers, genetic_id<Mode>);
        });
}

template <typename ModeDraws, typename JointDraws>
[[nodiscard]] auto make_marker_effects(
    std::type_identity<JointModeValues<ModeDraws, JointDraws>> /*draws*/,
    DrawReaders readers)
{
    auto mode_effects
        = make_marker_effects(std::type_identity<ModeDraws>{}, readers);
    auto joint_effects = make_marker_effects(
        std::type_identity<JointDraws>{}, readers, joint_genetic_id);
    return JointModeValues{std::move(mode_effects), std::move(joint_effects)};
}

template <typename Draws>
using genetic_marker_effects_t = decltype(make_marker_effects(
    std::type_identity<Draws>{},
    std::declval<DrawReaders>()));

template <GeneticModeSet Modes, typename... Effects, typename... Scales>
auto append_marker_columns(
    MarkerEffectTable& out,
    const ModeValues<Modes, Effects...>& effects,
    const ModeValues<Modes, Scales...>& scales) -> void
{
    effects.for_each(
        [&]<GeneticMode Mode>(const auto& mode_effects)
        {
            append_marker_columns(
                out, mode_effects, Mode, scales.template get<Mode>());
        });
}

// Per mode: coefficient columns, then that mode's PIP from the shared
// assignment; finally the any-effect PIP.
template <
    typename ModeEffects,
    typename JointEffects,
    GeneticModeSet Modes,
    typename... Scales>
auto append_marker_columns(
    MarkerEffectTable& out,
    const JointModeValues<ModeEffects, JointEffects>& effects,
    const ModeValues<Modes, Scales...>& scales) -> void
{
    effects.mode_values().for_each(
        [&]<GeneticMode Mode>(const auto& mode_effects)
        {
            append_marker_columns(
                out, mode_effects, Mode, scales.template get<Mode>());
            append_marker_columns(out, effects.joint(), Mode);
        });
    append_marker_columns(out, effects.joint());
}

template <typename Draws>
using genetic_diagnostics_t = decltype(make_diagnostics(
    std::type_identity<Draws>{},
    std::declval<DrawReaders>(),
    std::declval<double>()));

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_FACTORY_H_
