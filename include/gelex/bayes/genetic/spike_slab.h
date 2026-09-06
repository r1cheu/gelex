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
#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>

#include "gelex/bayes/basic_draw.h"
#include "gelex/bayes/basic_result.h"
#include "gelex/bayes/genetic/detail/draws_support.h"
#include "gelex/bayes/genetic/detail/pip_support.h"
#include "gelex/bayes/genetic/detail/prior_support.h"
#include "gelex/bayes/genetic/detail/result_support.h"
#include "gelex/bayes/genetic/detail/state_support.h"
#include "gelex/bayes/genetic/detail/summary_support.h"
#include "gelex/bayes/genetic/parameter.h"
#include "gelex/bayes/genetic/result.h"
#include "gelex/bayes/genetic/traits.h"
#include "gelex/bayes/genetic_policy.h"
#include "gelex/bayes/parameter.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/stats/dirichlet_log_kernel.h"
#include "gelex/bayes/variance/detail/calibration.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/detail/text_writer.h"
#include "gelex/namespace.h"

GELEX_NAMESPACE_BEGIN(gelex)

template <
    VarianceLayout Kind,
    MixtureWeightUpdate WeightUpdate = MixtureWeightUpdate::Enabled>
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
struct SpikeSlabPrior
{
    VarianceParameter variance;
    ProbabilityParameter<WeightUpdate> probability;
};

GELEX_NAMESPACE_BEGIN(detail)
template <
    GeneticMode Mode,
    VarianceLayout Kind,
    MixtureWeightUpdate WeightUpdate>
auto make_mode_prior(
    const SpikeSlabSpec<Kind, WeightUpdate>& spec,
    const MarkerVarianceCalibrator& calibrator)
    -> SpikeSlabPrior<Kind, WeightUpdate>
{
    return {
        .variance = calibrator.calibrate(Mode, spec.probability()),
        .probability = make_parameter<WeightUpdate>(
            spec.probability(), make_beta_prior(1.0, 1.0))};
}
GELEX_NAMESPACE_END(detail)

template <VarianceLayout Kind>
class SpikeSlabState
{
   public:
    SpikeSlabState(
        detail::marker_variance_state_t<Kind> variance,
        double probability,
        Eigen::Index num_markers,
        Eigen::Index num_individuals)
        : coefficients_(Eigen::VectorXd::Zero(num_markers)),
          assignments_(Eigen::VectorX<std::uint8_t>::Zero(num_markers)),
          class_counts_{static_cast<std::size_t>(num_markers), 0},
          fitted_values_(Eigen::VectorXd::Zero(num_individuals)),
          variance_(std::move(variance)),
          probability_(probability)
    {
    }

    auto coefficients() const -> const Eigen::VectorXd&
    {
        return coefficients_;
    }
    auto assignments() const -> const Eigen::VectorX<std::uint8_t>&
    {
        return assignments_;
    }
    auto class_counts() const -> const std::array<std::size_t, 2>&
    {
        return class_counts_;
    }
    auto fitted_values() const -> const Eigen::VectorXd&
    {
        return fitted_values_;
    }
    auto variance() const -> const detail::marker_variance_state_t<Kind>&
    {
        return variance_;
    }
    auto variance() -> detail::marker_variance_state_t<Kind>&
    {
        return variance_;
    }
    auto probability() const -> double { return probability_; }
    auto probability() -> double& { return probability_; }

    auto transition(Eigen::Index marker, double coefficient, bool active)
        -> void
    {
        const auto old_assignment = assignments_(marker);
        const auto new_assignment = static_cast<std::uint8_t>(active);
        if (old_assignment != new_assignment)
        {
            --class_counts_[old_assignment];
            ++class_counts_[new_assignment];
        }
        assignments_(marker) = new_assignment;
        coefficients_(marker) = active ? coefficient : 0.0;
    }

    auto transition(const Eigen::Ref<const Eigen::VectorXd>& delta) -> void
    {
        fitted_values_.noalias() += delta;
    }

   private:
    Eigen::VectorXd coefficients_;
    Eigen::VectorX<std::uint8_t> assignments_;
    std::array<std::size_t, 2> class_counts_;
    Eigen::VectorXd fitted_values_;
    detail::marker_variance_state_t<Kind> variance_;
    double probability_;
};

GELEX_NAMESPACE_BEGIN(detail)
template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
auto make_state(
    const SpikeSlabPrior<Kind, WeightUpdate>& prior,
    GeneticStateDimensions dimensions) -> SpikeSlabState<Kind>
{
    return {
        initial_marker_variance<Kind>(prior.variance, dimensions.marker_count),
        prior.probability.initial,
        dimensions.marker_count,
        dimensions.individual_count};
}
GELEX_NAMESPACE_END(detail)

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
struct SpikeSlabDraws
{
    detail::marker_variance_draw_t<Kind> variance;
    CategoryDraw<2> assignment;
    detail::weight_draw_t<WeightUpdate, ScalarDraw> probability;

    auto append(const SpikeSlabState<Kind>& state) -> void
    {
        variance.append(state.variance());
        assignment.append(state.assignments());
        probability.append(state.probability());
    }
};

GELEX_NAMESPACE_BEGIN(detail)
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
GELEX_NAMESPACE_END(detail)

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
struct SpikeSlabResult
{
    detail::marker_variance_result_t<Kind> variance;
    detail::weight_result_t<WeightUpdate, ScalarResult> probability;
};

GELEX_NAMESPACE_BEGIN(detail)
template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
auto make_result(const SpikeSlabDraws<Kind, WeightUpdate>& draws)
    -> SpikeSlabResult<Kind, WeightUpdate>
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
    const SpikeSlabResult<Kind, WeightUpdate>& result) -> void
{
    write_summary_rows(writer, result.variance);
    write_summary_rows(writer, result.probability);
}
GELEX_NAMESPACE_END(detail)

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_GENETIC_SPIKE_SLAB_H_
