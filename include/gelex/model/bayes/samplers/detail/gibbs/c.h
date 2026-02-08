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

#ifndef GELEX_MODEL_BAYES_SAMPLERS_DETAIL_GIBBS_C_H_
#define GELEX_MODEL_BAYES_SAMPLERS_DETAIL_GIBBS_C_H_

#include <random>

#include "../src/types/bayes_effects.h"
#include "gelex/model/bayes/samplers/detail/common_op.h"
#include "gelex/model/bayes/samplers/detail/gibbs/gibbs_concept.h"

namespace gelex::detail::Gibbs
{

template <typename EffectT, typename StateT>
    requires IsValidEffectStatePair<EffectT, StateT>
auto C(
    const EffectT& effect,
    StateT& state,
    bayes::ResidualState& residual,
    std::mt19937_64& rng) -> void
{
    auto& y_adj = residual.y_adj;
    const double residual_variance = residual.variance;

    const Eigen::VectorXd logpi = state.pi.prop.array().log();

    Eigen::VectorXd& coeffs = state.coeffs;
    auto& u = state.u;
    const double marker_variance = state.marker_variance(0);
    Eigen::VectorXi& tracker = state.tracker;

    const auto& X = bayes::get_matrix_ref(effect.X);
    const auto& cols_norm = effect.cols_norm;

    std::normal_distribution<double> normal{0, 1};
    std::uniform_real_distribution<double> uniform{0, 1};

    const double residual_over_marker_variance
        = residual_variance / marker_variance;

    double sum_square_coeffs{};

    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        if (effect.is_monomorphic(i))
        {
            continue;
        }

        const double old_i = coeffs(i);
        const auto& col = X.col(i);

        double rhs = blas_ddot(col, y_adj);
        if (old_i != 0.0)
        {
            rhs += cols_norm(i) * old_i;
        }

        auto [post_mean, post_stddev, log_like_kernel]
            = compute_posterior_params_core(
                rhs,
                cols_norm(i),
                residual_variance,
                residual_over_marker_variance);

        const double log_like_1_minus_0 = log_like_kernel + logpi(1) - logpi(0);

        const double prob_component_0
            = 1.0 / (1.0 + std::exp(log_like_1_minus_0));

        const int dist_index = (uniform(rng) < prob_component_0) ? 0 : 1;
        tracker(i) = dist_index;

        double new_i = 0.0;
        if (dist_index == 1)
        {
            new_i = (normal(rng) * post_stddev) + post_mean;
            update_residual_and_gebv(y_adj, u, col, old_i, new_i);
            sum_square_coeffs += new_i * new_i;
        }
        else if (old_i != 0.0)
        {
            update_residual_and_gebv(y_adj, u, col, old_i, 0.0);
        }
        coeffs(i) = new_i;
    }

    state.pi.count(1) = tracker.sum();
    state.pi.count(0) = static_cast<int>(coeffs.size() - state.pi.count(1));

    detail::ScaledInvChiSq chi_squared{effect.marker_variance_prior};
    chi_squared.compute(sum_square_coeffs, state.pi.count(1));
    state.marker_variance(0) = chi_squared(rng);

    state.variance = detail::var(state.u)(0);
}

}  // namespace gelex::detail::Gibbs

#endif  // GELEX_MODEL_BAYES_SAMPLERS_DETAIL_GIBBS_C_H_
