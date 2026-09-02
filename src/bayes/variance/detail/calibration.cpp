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

#include "gelex/bayes/variance/detail/calibration.h"

#include <cmath>
#include <fmt/format.h>

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

auto MarkerVarianceCalibrator::calibrate(
    GeneticMode mode,
    double initial_activity) const -> VarianceParameter
{
    const double projection_variance
        = model_->genetic().projection(mode).col_var().sum();
    if (projection_variance <= 0)
    {
        throw GelexException(
            fmt::format(
                "genetic column variance must sum to a positive value for mode "
                "{}",
                mode));
    }
    if (initial_activity <= 0)
    {
        throw GelexException(
            fmt::format(
                "initial marker activity must be positive, got {}",
                initial_activity));
    }

    const double target = model_->phenotype_variance() * budget_->share(mode)
                          / (initial_activity * projection_variance);
    return make_mean_calibrated_variance_parameter(target);
}

}  // namespace detail

}  // namespace gelex
