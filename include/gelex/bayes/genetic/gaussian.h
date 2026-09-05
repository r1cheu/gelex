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

#ifndef GELEX_BAYES_GENETIC_GAUSSIAN_H_
#define GELEX_BAYES_GENETIC_GAUSSIAN_H_

#include <Eigen/Core>

#include "gelex/bayes/basic_result.h"
#include "gelex/bayes/detail/genetic_spec.h"
#include "gelex/bayes/genetic/detail/draws_support.h"
#include "gelex/bayes/genetic/detail/result_support.h"
#include "gelex/bayes/genetic/detail/state_support.h"
#include "gelex/bayes/genetic/detail/summary_support.h"
#include "gelex/bayes/genetic/traits.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/parameter.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/variance/detail/calibration.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/detail/text_writer.h"

namespace gelex
{

template <VarianceLayout Kind>
struct GaussianFamily
{
};

template <VarianceLayout Kind>
struct GaussianPrior
{
    VarianceParameter variance;
};

template <VarianceLayout Kind>
struct GaussianState
{
    detail::marker_variance_state_t<Kind> variance{};
    Eigen::VectorXd fitted_values;  // total
};

template <VarianceLayout Kind>
[[nodiscard]] auto genetic_value(const GaussianState<Kind>& state)
    -> const Eigen::VectorXd&
{
    return state.fitted_values;
}

template <VarianceLayout Kind>
struct GaussianDraws
{
    detail::marker_variance_draw_t<Kind> variance;

    auto append(const GaussianState<Kind>& state) -> void
    {
        variance.append(state.variance);
    }
};

template <VarianceLayout Kind>
struct GaussianPosteriorResult
{
    detail::marker_variance_result_t<Kind> variance;
};

}  // namespace gelex

namespace gelex::detail
{

template <GeneticModeSet Modes, VarianceLayout Kind>
struct GeneticSpecFor<Modes, GaussianFamily<Kind>>
{
    using type = Gaussian;
};

template <GeneticModeSet Modes, VarianceLayout Kind>
auto make_prior(
    GaussianFamily<Kind> /*family*/,
    const Gaussian& /*genetic_spec*/,
    const MarkerVarianceCalibrator& calibrator)
{
    return generate_mode_values<Modes>(
        [&]<GeneticMode Mode>() -> GaussianPrior<Kind>
        { return {.variance = calibrator.calibrate(Mode, 1.0)}; });
}

template <VarianceLayout Kind>
auto make_state(
    const GaussianPrior<Kind>& prior,
    GeneticStateDimensions dimensions) -> GaussianState<Kind>
{
    return {
        .variance = initial_marker_variance<Kind>(
            prior.variance, dimensions.marker_count),
        .fitted_values = Eigen::VectorXd::Zero(dimensions.individual_count)};
}

template <VarianceLayout Kind>
[[nodiscard]] auto make_draws(
    const GaussianPrior<Kind>& /*prior*/,
    GeneticDrawsBuilder& builder) -> GaussianDraws<Kind>
{
    return {.variance = make_marker_variance_draw<Kind>(builder)};
}

template <VarianceLayout Kind>
auto make_result(const GaussianDraws<Kind>& draws)
    -> GaussianPosteriorResult<Kind>
{
    return {.variance = make_marker_variance_result<Kind>(draws.variance)};
}

template <VarianceLayout Kind>
[[nodiscard]] auto make_pip(const GaussianDraws<Kind>& /*draws*/)
    -> EmptyPosteriorResult
{
    return {};
}

template <VarianceLayout Kind>
auto write_family_summary_rows(
    TextWriter& writer,
    const GaussianPosteriorResult<Kind>& result) -> void
{
    write_summary_rows(writer, result.variance);
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_GAUSSIAN_H_
