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

#include "gelex/bayes/detail/calibration.h"

#include <cmath>
#include <fmt/format.h>

#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/variance_budget.h"
#include "gelex/exception.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{

namespace detail
{

namespace
{

constexpr double PRIOR_DEGREES_OF_FREEDOM = 4.0;

}  // namespace

auto MarkerVarianceCalibrator::calibrate(
    GeneticMode mode,
    double initial_activity) const -> VarianceParameter
{
    const double projection_variance = model_->genetic().col_var(mode).sum();
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
    if (!std::isfinite(target) || target <= 0)
    {
        throw GelexException(
            fmt::format(
                "calibrated marker variance must be finite and positive, got "
                "{}",
                target));
    }

    return VarianceParameter{
        .initial = target,
        .prior
        = {.degrees_of_freedom = PRIOR_DEGREES_OF_FREEDOM,
           .scale = (PRIOR_DEGREES_OF_FREEDOM - 2.0) / PRIOR_DEGREES_OF_FREEDOM
                    * target}};
}

}  // namespace detail

}  // namespace gelex
