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

#ifndef GELEX_MODEL_BAYES_SAMPLERS_DETAIL_GIBBS_A_H_
#define GELEX_MODEL_BAYES_SAMPLERS_DETAIL_GIBBS_A_H_

#include <random>

#include "../src/types/bayes_effects.h"
#include "gelex/model/bayes/samplers/detail/common_op.h"
#include "gelex/model/bayes/samplers/detail/gibbs/gibbs_concept.h"

namespace gelex::detail::Gibbs
{

template <typename EffectT, typename StateT>
    requires IsValidEffectStatePair<EffectT, StateT>
auto A(
    const EffectT& effect,
    StateT& state,
    bayes::ResidualState& residual,
    std::mt19937_64& rng) -> void
{
    auto& y_adj = residual.y_adj;
    const double residual_variance = residual.variance;

    Eigen::VectorXd& coeffs = state.coeffs;
    auto& u = state.u;
    Eigen::VectorXd& sigma = state.marker_variance;
    const auto& X = bayes::get_matrix_ref(effect.X);
    const auto& cols_norm = effect.cols_norm;

    detail::ScaledInvChiSq chi_squared{effect.marker_variance_prior};
    std::normal_distribution<double> normal{0, 1};

    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        if (effect.is_monomorphic(i))
        {
            continue;
        }
        const double old_i = coeffs(i);
        const auto& col = X.col(i);

        const double percision_kernel
            = 1 / (cols_norm(i) + residual_variance / sigma(i));

        // calculate the posterior mean and standard deviation
        const double rhs = blas_ddot(col, y_adj) + (cols_norm(i) * old_i);
        const double post_mean = rhs * percision_kernel;
        const double post_stddev = sqrt(residual_variance * percision_kernel);

        // sample a new coefficient
        const double new_i = (normal(rng) * post_stddev) + post_mean;
        coeffs(i) = new_i;

        chi_squared.compute(new_i * new_i);
        sigma(i) = chi_squared(rng);
        update_residual_and_gebv(y_adj, u, col, old_i, new_i);
    }
    state.variance = detail::var(state.u)(0);
}

}  // namespace gelex::detail::Gibbs

#endif  // GELEX_MODEL_BAYES_SAMPLERS_DETAIL_GIBBS_A_H_
