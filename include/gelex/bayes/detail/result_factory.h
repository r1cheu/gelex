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

#ifndef GELEX_BAYES_DETAIL_RESULT_FACTORY_H_
#define GELEX_BAYES_DETAIL_RESULT_FACTORY_H_

#include <Eigen/Core>
#include <cstddef>
#include <fmt/format.h>
#include <span>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "gelex/bayes/basic_draw.h"
#include "gelex/bayes/basic_result.h"
#include "gelex/bayes/genetic/draws.h"
#include "gelex/bayes/genetic/result.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

namespace gelex::detail
{

inline auto make_scalar_result(const ScalarDraw& draw) -> ScalarPosteriorResult
{
    return ScalarPosteriorResult{std::string{draw.identifier()}, draw.result()};
}

inline auto make_vector_result(
    const VectorDraw& draw,
    Eigen::Index expected_size) -> VectorPosteriorResult
{
    auto result
        = VectorPosteriorResult{std::string{draw.identifier()}, draw.result()};
    if (result.statistics().mean.size() != expected_size)
    {
        throw GelexException(
            fmt::format(
                "posterior '{}' has {} rows, expected {}",
                result.identifier(),
                result.statistics().mean.size(),
                expected_size));
    }
    return result;
}

inline auto make_coefficient_result(
    const VectorDraw& draw,
    std::span<const std::string> column_names,
    Eigen::Index expected_size) -> CoefficientPosteriorResult
{
    return CoefficientPosteriorResult{
        make_vector_result(draw, expected_size),
        std::vector<std::string>{column_names.begin(), column_names.end()}};
}

template <VarianceLayout Kind>
auto make_marker_variance_result(const marker_variance_draw_t<Kind>& draw)
    -> marker_variance_result_t<Kind>
{
    if constexpr (Kind == VarianceLayout::Pooled)
    {
        return make_scalar_result(draw);
    }
    else
    {
        static_cast<void>(draw);
        return {};
    }
}

template <UpdatePolicy Policy>
auto make_scalar_policy_result(const policy_draw_t<Policy, ScalarDraw>& draw)
    -> policy_result_t<Policy, ScalarPosteriorResult>
{
    if constexpr (Policy == UpdatePolicy::Sampled)
    {
        return make_scalar_result(draw);
    }
    else
    {
        static_cast<void>(draw);
        return {};
    }
}

template <UpdatePolicy Policy>
auto make_vector_policy_result(
    const policy_draw_t<Policy, VectorDraw>& draw,
    Eigen::Index expected_size)
    -> policy_result_t<Policy, VectorPosteriorResult>
{
    if constexpr (Policy == UpdatePolicy::Sampled)
    {
        return make_vector_result(draw, expected_size);
    }
    else
    {
        static_cast<void>(draw);
        static_cast<void>(expected_size);
        return {};
    }
}

template <VarianceLayout Kind>
auto make_family_result(const GaussianDraws<Kind>& draws)
    -> GaussianPosteriorResult<Kind>
{
    return {.variance = make_marker_variance_result<Kind>(draws.variance)};
}

inline auto make_family_result(
    const HalfNormalDraws<HalfNormalAsymmetry::Count>& draws)
    -> HalfNormalPosteriorResult<HalfNormalAsymmetry::Count>
{
    return {
        .variance = make_scalar_result(draws.variance),
        .positive_probability = make_scalar_result(draws.positive_probability)};
}

inline auto make_family_result(
    const HalfNormalDraws<HalfNormalAsymmetry::Magnitude>& draws)
    -> HalfNormalPosteriorResult<HalfNormalAsymmetry::Magnitude>
{
    return {.variances = make_vector_result(draws.variances, 2)};
}

template <VarianceLayout Kind, UpdatePolicy ProbabilityUpdate>
auto make_family_result(const SpikeSlabDraws<Kind, ProbabilityUpdate>& draws)
    -> SpikeSlabPosteriorResult<Kind, ProbabilityUpdate>
{
    return {
        .variance = make_marker_variance_result<Kind>(draws.variance),
        .probability
        = make_scalar_policy_result<ProbabilityUpdate>(draws.probability)};
}

template <std::size_t ClassCount, UpdatePolicy ProbabilitiesUpdate>
auto make_family_result(
    const ScaledMixtureDraws<ClassCount, ProbabilitiesUpdate>& draws)
    -> ScaledMixturePosteriorResult<ClassCount, ProbabilitiesUpdate>
{
    return {
        .variance = make_scalar_result(draws.variance),
        .probabilities = make_vector_policy_result<ProbabilitiesUpdate>(
            draws.probabilities, static_cast<Eigen::Index>(ClassCount)),
        .component_explained_variance = make_vector_result(
            draws.component_explained_variance,
            static_cast<Eigen::Index>(
                ScaledMixtureState<ClassCount>::component_count))};
}

template <std::size_t ClassCount, UpdatePolicy ProbabilitiesUpdate>
auto make_family_result(
    const JointSpikeSlabDraws<ClassCount, ProbabilitiesUpdate>& draws)
    -> JointSpikeSlabPosteriorResult<ClassCount, ProbabilitiesUpdate>
{
    return {
        .probabilities = make_vector_policy_result<ProbabilitiesUpdate>(
            draws.probabilities, static_cast<Eigen::Index>(ClassCount)),
        .component_explained_variance = make_vector_result(
            draws.component_explained_variance,
            static_cast<Eigen::Index>(
                JointSpikeSlabState<ClassCount>::component_count))};
}

template <GeneticModeSet Modes, typename... Priors, typename Draws>
auto make_genetic_result(
    std::type_identity<ModeValues<Modes, Priors...>> /*prior_type*/,
    const Draws& draws,
    Eigen::Index marker_count) -> genetic_result_t<ModeValues<Modes, Priors...>>
{
    return generate_mode_values<Modes>(
        [&]<GeneticMode Mode>()
        {
            const auto& mode_draws = draws.template get<Mode>();
            return GeneticModeResult{
                .coefficients
                = make_vector_result(mode_draws.coefficients, marker_count),
                .family_result = make_family_result(mode_draws.family_draws)};
        });
}

template <typename ModeValuesType, typename JointPrior, typename Draws>
auto make_genetic_result(
    std::type_identity<JointModeValues<ModeValuesType, JointPrior>>
    /*prior_type*/,
    const Draws& draws,
    Eigen::Index marker_count)
    -> genetic_result_t<JointModeValues<ModeValuesType, JointPrior>>
{
    auto mode_results = make_genetic_result(
        std::type_identity<ModeValuesType>{}, draws, marker_count);
    auto joint_result = make_family_result(draws.joint());
    return {std::move(mode_results), std::move(joint_result)};
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_DETAIL_RESULT_FACTORY_H_
