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

#ifndef GELEX_BAYES_GENETIC_FAMILY_H_
#define GELEX_BAYES_GENETIC_FAMILY_H_

// Lifts the per-method stage factories (make_prior, make_state, make_draws,
// make_result, make_pip, write_family_summary_rows) over ModeValues and
// JointModeValues so the top-level Bayes pipeline sees one genetic block.

#include <Eigen/Core>
#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <string>
#include <utility>

#include "gelex/bayes/basic_result.h"
#include "gelex/bayes/genetic/draws.h"
#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/joint_spike_slab.h"
#include "gelex/bayes/genetic/policy.h"
#include "gelex/bayes/genetic/result.h"
#include "gelex/bayes/genetic/scaled_mixture.h"
#include "gelex/bayes/genetic/spike_slab.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/variance/detail/calibration.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_writer.h"
#include "gelex/io/detail/text_writer.h"

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
        { return make_mode_prior<Mode>(spec, calibrator); });
}

// ---- state

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
        { return make_state(mode_prior, dimensions); });
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

// ---- draws

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
            GeneticDrawsBuilder builder{
                writer, fmt::format("genetic/{}", Mode), dimensions};
            return make_draws(mode_prior, builder);
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

// ---- pip

template <typename CoefficientDraws, typename ModeFamilyDraws>
[[nodiscard]] auto make_pip(
    const IndependentGeneticDraws<CoefficientDraws, ModeFamilyDraws>& draws)
{
    return generate_mode_values<CoefficientDraws::modes>(
        [&]<GeneticMode Mode>()
        { return make_pip(draws.template family<Mode>()); });
}

template <
    typename CoefficientDraws,
    typename ModeFamilyDraws,
    MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_pip(
    const JointGeneticDraws<
        CoefficientDraws,
        ModeFamilyDraws,
        JointSpikeSlabDraws<WeightUpdate>>& draws)
{
    const auto& assignment = draws.joint_family().assignment;
    auto mode_pip = generate_mode_values<GeneticMode::A | GeneticMode::D>(
        [&]<GeneticMode Mode>()
        {
            return MarkerPipResult{assignment.probability_of(
                [](std::size_t category)
                {
                    return JointSpikeSlabState::fitted_component_index<Mode>(
                               category)
                           != JointSpikeSlabState::no_component;
                })};
        });
    auto joint_pip
        = MarkerPipResult{assignment.probability_of(is_non_null_category)};
    return JointModeValues{std::move(mode_pip), std::move(joint_pip)};
}

// ---- result

template <typename CoefficientDraws, typename ModeFamilyDraws>
auto make_genetic_parameters(
    const IndependentGeneticDraws<CoefficientDraws, ModeFamilyDraws>& draws)
{
    return generate_mode_values<CoefficientDraws::modes>(
        [&]<GeneticMode Mode>()
        { return make_result(draws.template family<Mode>()); });
}

template <typename CoefficientDraws, typename ModeFamilyDraws, typename JointT>
auto make_genetic_parameters(
    const JointGeneticDraws<CoefficientDraws, ModeFamilyDraws, JointT>& draws)
{
    auto mode_results = generate_mode_values<CoefficientDraws::modes>(
        [&]<GeneticMode Mode>()
        { return make_result(draws.template family<Mode>()); });
    auto joint_result = make_result(draws.joint_family());
    return JointModeValues{std::move(mode_results), std::move(joint_result)};
}

template <typename GeneticDraws>
auto make_marker_effects(const BayesModel& model, const GeneticDraws& draws)
{
    constexpr auto modes = GeneticDraws::modes;
    const double phenotype_variance = model.phenotype_variance();
    auto pip = make_pip(draws);
    auto mode_results = generate_mode_values<modes>(
        [&]<GeneticMode Mode>()
        {
            const auto& coefficients
                = draws.coefficients().template get<Mode>();
            const auto& projection = model.genetic().projection(Mode);
            Eigen::VectorXd pve = projection.col_var().transpose().array()
                                  * coefficients.mean_square().array()
                                  / phenotype_variance;
            return MarkerEffectResult{
                make_result(coefficients),
                MarkerPveResult{std::move(pve)},
                std::move(pip.template get<Mode>())};
        });

    if constexpr (modes == (GeneticMode::A | GeneticMode::D))
    {
        Eigen::VectorXd joint_pve
            = mode_results.template get<GeneticMode::A>().pve().values()
              + mode_results.template get<GeneticMode::D>().pve().values();
        const auto covariance
            = model.genetic()
                  .projection(GeneticMode::A)
                  .col_covariance(model.genetic().projection(GeneticMode::D));
        joint_pve.array() += 2.0 * covariance.transpose().array()
                             * draws.coefficients().mean_product().array()
                             / phenotype_variance;

        auto joint_pip = [&]
        {
            if constexpr (requires { pip.joint(); })
            {
                return std::move(pip.joint());
            }
            else
            {
                return EmptyResult{};
            }
        }();
        return JointModeValues{
            std::move(mode_results),
            JointMarkerEffectResult{
                MarkerPveResult{std::move(joint_pve)}, std::move(joint_pip)}};
    }
    else
    {
        return mode_results;
    }
}

// ---- summary

template <GeneticModeSet Modes, typename... Results>
auto write_genetic_summary_rows(
    TextWriter& writer,
    const ModeValues<Modes, Results...>& result) -> void
{
    result.for_each([&]<GeneticMode Mode>(const auto& mode_result)
                    { write_family_summary_rows(writer, mode_result); });
}

template <typename ModeValuesType, typename JointResult>
auto write_genetic_summary_rows(
    TextWriter& writer,
    const JointModeValues<ModeValuesType, JointResult>& result) -> void
{
    write_genetic_summary_rows(writer, result.mode_values());
    write_family_summary_rows(writer, result.joint());
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_FAMILY_H_
