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
#include <cstddef>
#include <fmt/format.h>
#include <string_view>
#include <utility>

#include "gelex/bayes/genetic/draws.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/parameter.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/variance/detail/calibration.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_writer.h"
#include "gelex/namespace.h"

GELEX_NAMESPACE_BEGIN(gelex)

template <VarianceLayout Kind>
struct GaussianPrior
{
    VarianceParameter variance;
};

GELEX_NAMESPACE_BEGIN(detail)
template <GeneticMode Mode, VarianceLayout Kind>
auto make_prior(
    const GaussianSpec<Kind>& /*spec*/,
    const MarkerVarianceCalibrator& calibrator) -> GaussianPrior<Kind>
{
    return {.variance = calibrator.calibrate(Mode, 1.0)};
}
GELEX_NAMESPACE_END(detail)

template <VarianceLayout Kind>
class GaussianState
{
   public:
    GaussianState(
        detail::marker_variance_state_t<Kind> variance,
        GeneticDimensions dimensions)
        : coefficients_(
              Eigen::VectorXd::Zero(
                  static_cast<Eigen::Index>(dimensions.marker))),
          fitted_values_(
              Eigen::VectorXd::Zero(
                  static_cast<Eigen::Index>(dimensions.individual))),
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
auto make_state(const GaussianPrior<Kind>& prior, GeneticDimensions dimensions)
    -> GaussianState<Kind>
{
    return {
        initial_marker_variance<Kind>(
            prior.variance, static_cast<Eigen::Index>(dimensions.marker)),
        dimensions};
}
GELEX_NAMESPACE_END(detail)

template <VarianceLayout Kind>
class GaussianDraws
{
   public:
    using variance_writer_type = marker_variance_writer_t<Kind>;

    explicit GaussianDraws(
        variance_writer_type variances,
        PayloadWriter<float> coefficients)
        : variances_{std::move(variances)},
          coefficients_(std::move(coefficients))
    {
    }
    auto append(const GaussianState<Kind>& state) -> void
    {
        coefficients_.append(
            state.coefficients().template cast<float>().eval());
        if constexpr (Kind == VarianceLayout::Pooled)
        {
            variances_.append(state.variance());
        }
        else
        {
            variances_.append(state.variance().template cast<float>().eval());
        }
    }

   private:
    variance_writer_type variances_{};
    PayloadWriter<float> coefficients_;
};

GELEX_NAMESPACE_BEGIN(detail)
template <VarianceLayout Kind>
[[nodiscard]] auto make_draws(
    const GaussianPrior<Kind>& /*prior*/,
    BinaryWriter& writer,
    std::string_view prefix,
    std::size_t draw_count,
    GeneticDimensions dimensions) -> GaussianDraws<Kind>
{
    const std::size_t variance_size
        = (Kind == VarianceLayout::Pooled) ? 1 : dimensions.marker;

    auto variances = writer.reserve<marker_variance_dtype_t<Kind>>(
        fmt::format("{}/{}", prefix, variance_id),
        BinaryShape{variance_size, draw_count});
    auto coefficients = writer.reserve<float>(
        fmt::format("{}/{}", prefix, coefficients_id),
        BinaryShape{dimensions.marker, draw_count});

    return GaussianDraws<Kind>{std::move(variances), std::move(coefficients)};
}
GELEX_NAMESPACE_END(detail)
GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_GENETIC_GAUSSIAN_H_
