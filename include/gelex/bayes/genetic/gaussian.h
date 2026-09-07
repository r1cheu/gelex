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
#include "gelex/bayes/basic_result_io.h"
#include "gelex/bayes/genetic/draws.h"
#include "gelex/bayes/genetic/policy.h"
#include "gelex/bayes/genetic/result.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/parameter.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/variance/detail/calibration.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/detail/text_writer.h"
#include "gelex/namespace.h"

GELEX_NAMESPACE_BEGIN(gelex)

template <VarianceLayout Kind>
struct GaussianPrior
{
    VarianceParameter variance;
};

GELEX_NAMESPACE_BEGIN(detail)
template <GeneticModeSet Modes, VarianceLayout Kind>
auto make_prior(
    const GaussianSpec<Kind>& /*genetic_spec*/,
    const MarkerVarianceCalibrator& calibrator)
{
    return generate_mode_values<Modes>(
        [&]<GeneticMode Mode>() -> GaussianPrior<Kind>
        { return {.variance = calibrator.calibrate(Mode, 1.0)}; });
}
GELEX_NAMESPACE_END(detail)

template <VarianceLayout Kind>
class GaussianState
{
   public:
    GaussianState(
        detail::marker_variance_state_t<Kind> variance,
        Eigen::Index num_markers,
        Eigen::Index num_individuals)
        : coefficients_(Eigen::VectorXd::Zero(num_markers)),
          fitted_values_(Eigen::VectorXd::Zero(num_individuals)),
          variance_(variance)
    {
    }

    auto coefficients() const -> const Eigen::VectorXd&
    {
        return coefficients_;
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

    auto transition(Eigen::Index marker_index, double coefficient) -> void
    {
        coefficients_(marker_index) = coefficient;
    }
    auto transition(const Eigen::Ref<const Eigen::VectorXd>& delta) -> void
    {
        fitted_values_.noalias() += delta;
    }

   private:
    Eigen::VectorXd coefficients_;
    Eigen::VectorXd fitted_values_;
    detail::marker_variance_state_t<Kind> variance_;
};

GELEX_NAMESPACE_BEGIN(detail)
template <VarianceLayout Kind>
auto make_state(
    const GaussianPrior<Kind>& prior,
    GeneticStateDimensions dimensions) -> GaussianState<Kind>
{
    return {
        initial_marker_variance<Kind>(prior.variance, dimensions.marker_count),
        dimensions.marker_count,
        dimensions.individual_count};
}
GELEX_NAMESPACE_END(detail)

template <VarianceLayout Kind>
struct GaussianDraws
{
    detail::marker_variance_draw_t<Kind> variance;

    auto append(const GaussianState<Kind>& state) -> void
    {
        variance.append(state.variance());
    }
};

GELEX_NAMESPACE_BEGIN(detail)
template <VarianceLayout Kind>
[[nodiscard]] auto make_draws(
    const GaussianPrior<Kind>& /*prior*/,
    GeneticDrawsBuilder& builder) -> GaussianDraws<Kind>
{
    return {.variance = make_marker_variance_draw<Kind>(builder)};
}
GELEX_NAMESPACE_END(detail)

template <VarianceLayout Kind>
struct GaussianResult
{
    detail::marker_variance_result_t<Kind> variance;
};

GELEX_NAMESPACE_BEGIN(detail)
template <VarianceLayout Kind>
auto make_result(const GaussianDraws<Kind>& draws) -> GaussianResult<Kind>
{
    return {.variance = make_marker_variance_result<Kind>(draws.variance)};
}

template <VarianceLayout Kind>
[[nodiscard]] auto make_pip(const GaussianDraws<Kind>& /*draws*/) -> EmptyResult
{
    return {};
}

template <VarianceLayout Kind>
auto write_family_summary_rows(
    TextWriter& writer,
    const GaussianResult<Kind>& result) -> void
{
    write_summary_rows(writer, result.variance);
}
GELEX_NAMESPACE_END(detail)

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_GENETIC_GAUSSIAN_H_
