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

#ifndef GELEX_ALGO_INFER_VI_KERNELS_RR_H_
#define GELEX_ALGO_INFER_VI_KERNELS_RR_H_

#include <variant>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/state.h"
#include "gelex/algo/infer/vi/kernels/concept.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/legacy_method.h"

namespace gelex::vi
{

// CAVI ridge: closed-form Gaussian posterior per marker, ScaledInvChiSq
// expectation for the shared marker variance.
// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class RRKernel
{
   public:
    RRKernel(
        const bayes::OldGeneticPrior& prior,
        const bayes::GeneticEffect& /*effect*/,
        bayes::vi::GeneticState& state)
        : state_(state),
          chi_squared_(
              make_chi_squared_params(std::get<bayes::GeneticSpec>(prior.spec)))
    {
    }

    auto prepare() -> void
    {
        inv_marker_variance_ = 1.0 / state_.marker_variance(0);
        n_used_ = 0;
    }

    auto update(
        Eigen::Index /*marker_index*/,
        double xtx_diag,
        double rhs,
        double residual_variance) -> CaviUpdate
    {
        const double v = xtx_diag + (residual_variance * inv_marker_variance_);
        const double inv_v = 1.0 / v;
        ++n_used_;
        return {.mean = rhs * inv_v, .sigma2 = residual_variance * inv_v};
    }

    auto commit() -> void
    {
        const double sum_sq = state_.coeffs.squaredNorm() + state_.sigma2.sum();
        chi_squared_.compute(sum_sq, n_used_);
        state_.marker_variance(0) = chi_squared_.expected_value();
    }

   private:
    static auto make_chi_squared_params(const bayes::GeneticSpec& spec)
        -> stats::detail::ScaledInvChiSqParams
    {
        return {
            .nu = spec.variance.prior.degrees_of_freedom,
            .s2 = spec.variance.prior.scale};
    }

    bayes::vi::GeneticState& state_;
    stats::detail::ScaledInvChiSq chi_squared_;
    double inv_marker_variance_{};
    Eigen::Index n_used_{0};
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::vi

#endif  // GELEX_ALGO_INFER_VI_KERNELS_RR_H_
