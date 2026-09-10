// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_GAUSSIAN_H_
#define GELEX_BAYES_GENETIC_GAUSSIAN_H_

#include <Eigen/Core>
#include <cassert>
#include <cstddef>
#include <fmt/format.h>
#include <string_view>
#include <utility>

#include "gelex/bayes/genetic/detail/marker_variance.h"
#include "gelex/bayes/genetic/draw_traits.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/parameter.h"
#include "gelex/bayes/serialization_ids.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/variance/calibration.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/dense_writer.h"
#include "gelex/namespace.h"

GELEX_NAMESPACE_BEGIN(gelex)

template <VarianceLayout Kind>
struct GaussianPrior
{
    VarianceParameter variance;
};

template <GeneticMode Mode, VarianceLayout Kind>
auto make_prior(
    const GaussianSpec<Kind>& /*spec*/,
    const MarkerVarianceCalibrator& calibrator) -> GaussianPrior<Kind>
{
    return {.variance = calibrator.calibrate(Mode, 1.0)};
}

template <VarianceLayout Kind>
class GaussianState
{
   public:
    using variance_type = detail::marker_variance_state_t<Kind>;
    GaussianState(variance_type variance, GeneticDimensions dimensions)
        : coefficients_(
              Eigen::VectorXd::Zero(
                  static_cast<Eigen::Index>(dimensions.marker))),
          variance_(variance)
    {
    }

    auto coefficients() const -> const Eigen::VectorXd&
    {
        return coefficients_;
    }
    auto variance() const -> const variance_type& { return variance_; }
    auto variance() -> variance_type& { return variance_; }

    auto transition(Eigen::Index marker, double coefficient) -> void
    {
        assert(marker >= 0 && marker < coefficients_.size());
        coefficients_(marker) = coefficient;
    }

   private:
    Eigen::VectorXd coefficients_;
    variance_type variance_;
};

template <VarianceLayout Kind>
auto make_state(const GaussianPrior<Kind>& prior, GeneticDimensions dimensions)
    -> GaussianState<Kind>
{
    return {
        detail::initial_marker_variance<Kind>(
            prior.variance, static_cast<Eigen::Index>(dimensions.marker)),
        dimensions};
}

template <VarianceLayout Kind, CoefficientLayout Layout>
class GaussianDraws
{
   public:
    using variance_writer_type = marker_variance_writer_t<Kind>;
    using coefficients_writer_type = coefficients_writer_t<Layout>;

    explicit GaussianDraws(
        variance_writer_type variances,
        coefficients_writer_type coefficients)
        : variances_{std::move(variances)},
          coefficients_(std::move(coefficients))
    {
    }
    auto operator<<(const GaussianState<Kind>& state) -> GaussianDraws&
    {
        coefficients_ << state.coefficients();
        if constexpr (Kind == VarianceLayout::Pooled)
        {
            variances_ << state.variance();
        }
        return *this;
    }

   private:
    [[no_unique_address]] variance_writer_type variances_;
    coefficients_writer_type coefficients_;
};

template <
    CoefficientLayout Layout = CoefficientLayout::Dense,
    VarianceLayout Kind>
[[nodiscard]] auto make_draws(
    const GaussianState<Kind>& state,
    DrawWriters writers,
    std::string_view prefix,
    std::size_t draw_count) -> GaussianDraws<Kind, Layout>
{
    const auto marker_count
        = static_cast<std::size_t>(state.coefficients().size());

    auto variances = reserve_marker_variance<Kind>(
        writers, fmt::format("{}/{}", prefix, variance_id), draw_count);
    auto coefficients = reserve_coefficients<Layout>(
        writers,
        fmt::format("{}/{}", prefix, coefficients_id),
        BinaryShape{marker_count, draw_count});

    return GaussianDraws<Kind, Layout>{
        std::move(variances), std::move(coefficients)};
}

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_GENETIC_GAUSSIAN_H_
