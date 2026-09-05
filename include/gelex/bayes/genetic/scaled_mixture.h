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

#ifndef GELEX_BAYES_GENETIC_SCALED_MIXTURE_H_
#define GELEX_BAYES_GENETIC_SCALED_MIXTURE_H_

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <ranges>

#include "gelex/bayes/basic_draw.h"
#include "gelex/bayes/basic_result.h"
#include "gelex/bayes/detail/genetic_spec.h"
#include "gelex/bayes/genetic/detail/draws_support.h"
#include "gelex/bayes/genetic/detail/pip_support.h"
#include "gelex/bayes/genetic/detail/prior_support.h"
#include "gelex/bayes/genetic/detail/result_support.h"
#include "gelex/bayes/genetic/detail/state_support.h"
#include "gelex/bayes/genetic/detail/summary_support.h"
#include "gelex/bayes/genetic/parameter.h"
#include "gelex/bayes/genetic/result.h"
#include "gelex/bayes/genetic/traits.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/parameter.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/stats/dirichlet_log_kernel.h"
#include "gelex/bayes/variance/detail/calibration.h"
#include "gelex/genetic_mode.h"
#include "gelex/infra/var.h"
#include "gelex/io/detail/text_writer.h"

namespace gelex
{

template <MixtureWeightUpdate WeightUpdate = MixtureWeightUpdate::Enabled>
struct ScaledMixtureFamily
{
};

template <
    std::size_t ClassCount,
    MixtureWeightUpdate WeightUpdate = MixtureWeightUpdate::Enabled>
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
struct ScaledMixturePrior
{
    static constexpr std::size_t class_count = ClassCount;

    VarianceParameter variance;
    SimplexParameter<class_count, WeightUpdate> probabilities;
    std::array<double, class_count> scales{};
};

template <std::size_t ClassCount>
    requires(
        ClassCount > 1
        && ClassCount <= static_cast<std::size_t>(
                             std::numeric_limits<std::uint8_t>::max())
                             + 1)
struct ScaledMixtureState
{
    static constexpr std::size_t class_count = ClassCount;
    static constexpr std::size_t component_count = class_count - 1;

    double variance{};
    Eigen::VectorX<std::uint8_t> assignment;
    std::array<double, class_count> probabilities{};
    Eigen::Matrix<double, Eigen::Dynamic, static_cast<int>(component_count)>
        fitted_values;
};

template <std::size_t ClassCount>
[[nodiscard]] auto genetic_value(const ScaledMixtureState<ClassCount>& state)
{
    return state.fitted_values.rowwise().sum();
}

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
struct ScaledMixtureDraws
{
    ScalarDraw variance;
    CategoryDraw<ClassCount> assignment;
    detail::weight_draw_t<WeightUpdate, VectorDraw> probabilities;
    VectorDraw component_explained_variance;

    auto append(const ScaledMixtureState<ClassCount>& state) -> void
    {
        variance.append(state.variance);
        assignment.append(state.assignment);
        probabilities.append(state.probabilities);
        component_explained_variance.append(
            matvar<0>(state.fitted_values, VarNormType::Population));
    }
};

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
struct ScaledMixtureResult
{
    ScalarResult variance;
    detail::weight_result_t<WeightUpdate, VectorResult> probabilities;
    VectorResult component_explained_variance;
};

}  // namespace gelex

namespace gelex::detail
{

template <GeneticModeSet Modes, MixtureWeightUpdate WeightUpdate>
struct GeneticSpecFor<Modes, ScaledMixtureFamily<WeightUpdate>>
{
    using type = HomogeneousModeValues<Modes, ScaledMixture>;
};

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

template <GeneticModeSet Modes, MixtureWeightUpdate WeightUpdate>
auto make_prior(
    ScaledMixtureFamily<WeightUpdate> /*family*/,
    const genetic_spec_t<Modes, ScaledMixtureFamily<WeightUpdate>>&
        genetic_spec,
    const MarkerVarianceCalibrator& calibrator)
{
    return transform_mode_values(
        genetic_spec,
        [&]<GeneticMode Mode>(const ScaledMixture& spec)
            -> ScaledMixturePrior<ScaledMixture::class_count, WeightUpdate>
        {
            return {
                .variance = calibrator.calibrate(Mode, initial_activity(spec)),
                .probabilities = make_parameter<WeightUpdate>(
                    spec.probabilities(),
                    make_uniform_dirichlet_prior<ScaledMixture::class_count>()),
                .scales = spec.scales()};
        });
}

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
auto make_state(
    const ScaledMixturePrior<ClassCount, WeightUpdate>& prior,
    GeneticStateDimensions dimensions) -> ScaledMixtureState<ClassCount>
{
    auto state = ScaledMixtureState<ClassCount>{
        .variance = prior.variance.initial,
        .assignment
        = Eigen::VectorX<std::uint8_t>::Zero(dimensions.marker_count),
        .probabilities = prior.probabilities.initial,
    };
    state.fitted_values.setZero(
        dimensions.individual_count,
        static_cast<Eigen::Index>(state.component_count));
    return state;
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

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
auto make_result(const ScaledMixtureDraws<ClassCount, WeightUpdate>& draws)
    -> ScaledMixtureResult<ClassCount, WeightUpdate>
{
    return {
        .variance = make_result(draws.variance),
        .probabilities = make_result(draws.probabilities),
        .component_explained_variance
        = make_result(draws.component_explained_variance)};
}

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_pip(
    const ScaledMixtureDraws<ClassCount, WeightUpdate>& draws)
    -> MarkerPipResult
{
    return MarkerPipResult{
        draws.assignment.probability_of(is_non_null_category)};
}

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
auto write_family_summary_rows(
    TextWriter& writer,
    const ScaledMixtureResult<ClassCount, WeightUpdate>& result) -> void
{
    write_summary_rows(writer, result.variance);
    write_summary_rows(writer, result.probabilities);
    write_summary_rows(writer, result.component_explained_variance);
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_SCALED_MIXTURE_H_
