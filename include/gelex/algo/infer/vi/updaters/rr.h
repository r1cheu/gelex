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

#ifndef GELEX_ALGO_INFER_VI_UPDATERS_RR_H_
#define GELEX_ALGO_INFER_VI_UPDATERS_RR_H_

#include "gelex/algo/infer/detail/marker_op.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/infra/stats/descriptive.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/states.h"
#include "gelex/model/bayes/vi/states.h"

namespace gelex::vi::detail
{

inline auto RR(
    const bayes::GeneticEffect& effect,
    const bayes::GeneticPrior& prior,
    bayes::vi::GeneticState& state,
    bayes::ResidualState& residual) -> void
{
    auto& y_adj = residual.y_adj;
    const double residual_variance = residual.variance;

    Eigen::VectorXd& coeffs = state.coeffs;
    Eigen::VectorXd& sigma2 = state.sigma2;
    const double old_marker_variance = state.marker_variance(0);
    Eigen::VectorXd& u = state.u;
    const auto& X = bayes::get_matrix_ref(effect.X);
    const auto& XtX_diag = effect.XtX_diag;

    const double residual_over_var = residual_variance / old_marker_variance;

    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        if (effect.is_monomorphic(i))
        {
            continue;
        }

        const double old_i = coeffs(i);
        const auto col = X.col(i);
        const double v = XtX_diag(i) + residual_over_var;
        const double inv_v = 1.0 / v;

        const double rhs
            = gelex::detail::blas_ddot(col, y_adj) + (XtX_diag(i) * old_i);
        const double post_mean = rhs * inv_v;
        const double post_var = residual_variance * inv_v;

        coeffs(i) = post_mean;
        sigma2(i) = post_var;
        gelex::detail::apply_marker_update(
            y_adj, u, {}, col, {.old_value = old_i, .new_value = post_mean});
    }
    state.variance = gelex::detail::var(state.u)(0);

    const auto& marker_prior = std::get<bayes::ContinuousPrior>(prior.marker);
    gelex::detail::ScaledInvChiSq chi_squared{marker_prior.variance.param};
    // E[beta^2] = mu^2 + sigma2, unlike MCMC which uses sampled beta^2
    const double sum_sq = (coeffs.array().square() + sigma2.array()).sum();
    chi_squared.compute(sum_sq, coeffs.size() - effect.num_mono());
    state.marker_variance(0) = chi_squared.expected_value();
}

}  // namespace gelex::vi::detail

#endif  // GELEX_ALGO_INFER_VI_UPDATERS_RR_H_
