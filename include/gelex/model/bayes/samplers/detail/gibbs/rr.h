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

#ifndef GELEX_MODEL_BAYES_SAMPLERS_DETAIL_GIBBS_RR_H_
#define GELEX_MODEL_BAYES_SAMPLERS_DETAIL_GIBBS_RR_H_

#include <random>

#include "../src/types/bayes_effects.h"
#include "gelex/model/bayes/samplers/detail/common_op.h"
#include "gelex/model/bayes/samplers/detail/gibbs/gibbs_concept.h"

namespace gelex::detail::Gibbs
{

template <typename EffectT, typename StateT>
    requires IsValidEffectStatePair<EffectT, StateT>
auto RR(
    const EffectT& effect,
    StateT& state,
    bayes::ResidualState& residual,
    std::mt19937_64& rng) -> void
{
    auto& y_adj = residual.y_adj;
    const double residual_variance = residual.variance;

    Eigen::VectorXd& coeff = state.coeffs;
    const double old_marker_variance = state.marker_variance(0);
    Eigen::VectorXd& u = state.u;
    const auto& X = bayes::get_matrix_ref(effect.X);
    const auto& cols_norm = effect.cols_norm;

    const double residual_over_var = residual_variance / old_marker_variance;
    const double sqrt_residual_variance = std::sqrt(residual_variance);

    std::normal_distribution<double> normal{0, 1};

    for (Eigen::Index i = 0; i < coeff.size(); ++i)
    {
        if (effect.is_monomorphic(i))
        {
            continue;
        }

        const double old_i = coeff(i);
        const auto col = X.col(i);
        const double v = cols_norm(i) + residual_over_var;
        const double inv_v = 1.0 / v;

        const double rhs = blas_ddot(col, y_adj) + (cols_norm(i) * old_i);
        const double post_mean = rhs * inv_v;
        const double post_stddev = sqrt_residual_variance * std::sqrt(inv_v);

        const double new_i = (normal(rng) * post_stddev) + post_mean;
        coeff(i) = new_i;
        update_residual_and_gebv(y_adj, u, col, old_i, new_i);
    }
    state.variance = detail::var(state.u)(0);

    detail::ScaledInvChiSq chi_squared{effect.marker_variance_prior};
    chi_squared.compute(coeff.squaredNorm(), coeff.size() - effect.num_mono());
    state.marker_variance(0) = chi_squared(rng);
}

}  // namespace gelex::detail::Gibbs

#endif  // GELEX_MODEL_BAYES_SAMPLERS_DETAIL_GIBBS_RR_H_
