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

#ifndef GELEX_ALGO_INFER_MCMC_SWEEPS_FIXED_COEFFICIENT_H_
#define GELEX_ALGO_INFER_MCMC_SWEEPS_FIXED_COEFFICIENT_H_

#include <random>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/invariant.h"
#include "gelex/model/bayes/state.h"
#include "gelex/types/fixed_designs.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class FixedCoefficientSweep
{
   public:
    FixedCoefficientSweep(
        const FixedDesign& design,
        bayes::FixedState& state,
        bayes::ResidualState& residual)
        : design_(design), state_(state), residual_(residual)
    {
    }

    template <typename Kernel>
    auto run(Kernel& kernel, std::mt19937_64& rng) -> void
    {
        const auto& X = design_.X;
        const auto& XtX_diag = design_.XtX_diag;
        auto& coeffs = state_.coeffs;
        const double residual_variance = residual_.variance;

        kernel.prepare();
        for (Eigen::Index i = 0; i < coeffs.size(); ++i)
        {
            const auto col = X.col(i);
            const double old_i = coeffs(i);
            ResidualAdjustmentGuard guard{col, coeffs(i), residual_};
            const double rhs = col.dot(residual_.y_adj) + (XtX_diag(i) * old_i);
            coeffs(i) = kernel.sample(XtX_diag(i), rhs, residual_variance, rng);
        }
    }

   private:
    const FixedDesign& design_;
    bayes::FixedState& state_;
    bayes::ResidualState& residual_;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_SWEEPS_FIXED_COEFFICIENT_H_
