// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/bayes/variance/calibration.h"

#include <array>
#include <cmath>
#include <fmt/format.h>
#include <utility>

#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/genotype/projection.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/stats/scaled_inv_chi2_log_kernel.h"
#include "gelex/bayes/variance/budget.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

namespace gelex
{

namespace detail
{

namespace
{

constexpr double prior_degrees_of_freedom = 4.0;

}  // namespace

auto make_mean_calibrated_variance_parameter(double target) -> VarianceParameter
{
    if (!std::isfinite(target) || target <= 0.0)
    {
        throw GelexException(
            fmt::format(
                "variance target must be finite and positive, got {}", target));
    }

    return VarianceParameter{
        .initial = target,
        .prior = make_scaled_inv_chi2_prior(
            prior_degrees_of_freedom,
            (prior_degrees_of_freedom - 2.0) / prior_degrees_of_freedom
                * target)};
}

}  // namespace detail

MarkerVarianceCalibrator::MarkerVarianceCalibrator(
    std::array<double, all_genetic_modes.size()> base_variances)
    : base_variances_{base_variances}
{
    for (const double variance : base_variances_)
    {
        if (!std::isfinite(variance) || variance < 0.0)
        {
            throw GelexException(
                "base marker variance must be finite and non-negative");
        }
    }
}

auto make_marker_variance_calibrator(
    const BayesModel& model,
    const VarianceBudget& budget) -> MarkerVarianceCalibrator
{
    std::array<double, all_genetic_modes.size()> base_variances{};
    for (const auto mode : model.genetic().each_mode())
    {
        const double projection_variance
            = model.genetic().projection(mode).col_var().sum();
        if (!std::isfinite(projection_variance) || projection_variance <= 0.0)
        {
            throw GelexException(
                fmt::format(
                    "genetic column variance must sum to a finite positive "
                    "value for mode {}",
                    mode));
        }
        base_variances[std::to_underlying(mode)] = model.phenotype_variance()
                                                   * budget.genetic(mode)
                                                   / projection_variance;
    }
    return MarkerVarianceCalibrator{base_variances};
}

auto MarkerVarianceCalibrator::calibrate(
    GeneticMode mode,
    double initial_activity) const -> VarianceParameter
{
    if (!std::isfinite(initial_activity) || initial_activity <= 0.0)
    {
        throw GelexException(
            fmt::format(
                "initial marker activity must be finite and positive, got {}",
                initial_activity));
    }
    const double target
        = base_variances_[std::to_underlying(mode)] / initial_activity;
    return detail::make_mean_calibrated_variance_parameter(target);
}

}  // namespace gelex
