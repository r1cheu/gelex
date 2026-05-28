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

#include "gelex/algo/infer/mcmc/steps/fixed_coefficient.h"

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/invariant.h"
#include "gelex/infra/stats/conjugate_prior.h"

namespace gelex::mcmc
{

FixedCoefficientStep::FixedCoefficientStep(
    const FixedDesign& design,
    bayes::FixedState& state,
    bayes::ResidualState& residual,
    std::mt19937_64& rng)
    : design_(design), state_(state), residual_(residual), rng_(rng)
{
}

auto FixedCoefficientStep::step() -> void
{
    const double residual_variance = residual_.variance;
    auto& coeffs = state_.coeffs;

    normal_.reset();
    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        const auto column = design_.X.col(i);
        const double old_value = coeffs(i);
        ResidualAdjustmentGuard guard{column, coeffs(i), residual_};
        const double xtx_diag_i = design_.XtX_diag(i);
        const double rhs
            = column.dot(residual_.y_adj) + (xtx_diag_i * old_value);

        coeffs(i) = normal_.draw(
            {.mean = rhs / xtx_diag_i, .var = residual_variance / xtx_diag_i},
            rng_);
    }
}

}  // namespace gelex::mcmc
