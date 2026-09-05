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

#ifndef GELEX_BAYES_GENETIC_SPIKE_SLAB_H_
#define GELEX_BAYES_GENETIC_SPIKE_SLAB_H_

#include <Eigen/Core>
#include <cstdint>

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
#include "gelex/io/detail/text_writer.h"

namespace gelex
{

template <
    VarianceLayout Kind,
    MixtureWeightUpdate WeightUpdate = MixtureWeightUpdate::Enabled>
struct SpikeSlabFamily
{
};

template <
    VarianceLayout Kind,
    MixtureWeightUpdate WeightUpdate = MixtureWeightUpdate::Enabled>
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
struct SpikeSlabPrior
{
    VarianceParameter variance;
    ProbabilityParameter<WeightUpdate> probability;
};

template <VarianceLayout Kind>
struct SpikeSlabState
{
    detail::marker_variance_state_t<Kind> variance{};
    Eigen::VectorX<std::uint8_t> assignment;
    double probability{};
    Eigen::VectorXd fitted_values;  // total
};

template <VarianceLayout Kind>
[[nodiscard]] auto genetic_value(const SpikeSlabState<Kind>& state)
    -> const Eigen::VectorXd&
{
    return state.fitted_values;
}

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
struct SpikeSlabDraws
{
    detail::marker_variance_draw_t<Kind> variance;
    CategoryDraw<2> assignment;
    detail::weight_draw_t<WeightUpdate, ScalarDraw> probability;

    auto append(const SpikeSlabState<Kind>& state) -> void
    {
        variance.append(state.variance);
        assignment.append(state.assignment);
        probability.append(state.probability);
    }
};

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
struct SpikeSlabPosteriorResult
{
    detail::marker_variance_result_t<Kind> variance;
    detail::weight_result_t<WeightUpdate, ScalarPosteriorResult> probability;
};

}  // namespace gelex

namespace gelex::detail
{

template <
    GeneticModeSet Modes,
    VarianceLayout Kind,
    MixtureWeightUpdate WeightUpdate>
struct GeneticSpecFor<Modes, SpikeSlabFamily<Kind, WeightUpdate>>
{
    using type = HomogeneousModeValues<Modes, SpikeSlab>;
};

template <
    GeneticModeSet Modes,
    VarianceLayout Kind,
    MixtureWeightUpdate WeightUpdate>
auto make_prior(
    SpikeSlabFamily<Kind, WeightUpdate> /*family*/,
    const genetic_spec_t<Modes, SpikeSlabFamily<Kind, WeightUpdate>>&
        genetic_spec,
    const MarkerVarianceCalibrator& calibrator)
{
    return transform_mode_values(
        genetic_spec,
        [&]<GeneticMode Mode>(
            const SpikeSlab& spec) -> SpikeSlabPrior<Kind, WeightUpdate>
        {
            return {
                .variance = calibrator.calibrate(Mode, spec.probability()),
                .probability = make_parameter<WeightUpdate>(
                    spec.probability(), make_beta_prior(1.0, 1.0))};
        });
}

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
auto make_state(
    const SpikeSlabPrior<Kind, WeightUpdate>& prior,
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

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
auto make_result(const SpikeSlabDraws<Kind, WeightUpdate>& draws)
    -> SpikeSlabPosteriorResult<Kind, WeightUpdate>
{
    return {
        .variance = make_marker_variance_result<Kind>(draws.variance),
        .probability = make_result(draws.probability)};
}

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_pip(const SpikeSlabDraws<Kind, WeightUpdate>& draws)
    -> MarkerPipResult
{
    return MarkerPipResult{
        draws.assignment.probability_of(is_non_null_category)};
}

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
auto write_family_summary_rows(
    TextWriter& writer,
    const SpikeSlabPosteriorResult<Kind, WeightUpdate>& result) -> void
{
    write_summary_rows(writer, result.variance);
    write_summary_rows(writer, result.probability);
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_SPIKE_SLAB_H_
