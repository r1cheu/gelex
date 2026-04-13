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

#ifndef GELEX_ALGO_INFER_MCMC_SAMPLERS_GIBBS_B_H_
#define GELEX_ALGO_INFER_MCMC_SAMPLERS_GIBBS_B_H_

#include <random>

#include "gelex/algo/infer/mcmc/samplers/common_op.h"
#include "gelex/algo/infer/mcmc/samplers/gibbs/gibbs_concept.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/states.h"

namespace gelex::detail::Gibbs
{

template <typename EffectT, typename StateT>
    requires IsValidEffectStatePair<EffectT, StateT>
auto B(
    const EffectT& effect,
    const bayes::GeneticPrior& prior,
    StateT& state,
    bayes::ResidualState& residual,
    std::mt19937_64& rng) -> void
{
    auto& y_adj = residual.y_adj;
    const double residual_variance = residual.variance;

    auto& marker_assignment = std::get<bayes::Assignment>(*state.group);
    const double logpi_0 = std::log(marker_assignment.proportion(0));
    const double logpi_1 = std::log(marker_assignment.proportion(1));

    Eigen::VectorXd& coeffs = state.coeffs;
    auto& u = state.u;
    Eigen::VectorXd& marker_variance = state.marker_variance;
    auto& tracker = marker_assignment.tracker;

    const auto& X = bayes::get_matrix_ref(effect.X);
    const auto& cols_squared_norm = effect.cols_squared_norm;

    std::normal_distribution<double> normal{0, 1};
    std::uniform_real_distribution<double> uniform{0, 1};
    const auto& marker_prior = std::get<bayes::SpikePrior>(prior.marker);
    detail::ScaledInvChiSq chi_squared{marker_prior.variance.param};

    int count_1 = 0;
    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        if (effect.is_monomorphic(i))
        {
            count_1 += tracker(i);
            continue;
        }

        const double old_i = coeffs(i);
        const auto& col = X.col(i);
        const double variance_i = marker_variance(i);

        double rhs = blas_ddot(col, y_adj);
        if (old_i != 0.0)
        {
            rhs += cols_squared_norm(i) * old_i;
        }

        auto [post_mean, post_stddev, log_like_kernel]
            = compute_posterior_params(
                rhs, variance_i, cols_squared_norm(i), residual_variance);

        const double log_like_1_minus_0 = log_like_kernel + logpi_1 - logpi_0;

        const double prob_component_0
            = 1.0 / (1.0 + std::exp(log_like_1_minus_0));

        const int8_t dist_index = (uniform(rng) < prob_component_0) ? 0 : 1;
        tracker(i) = dist_index;
        count_1 += dist_index;

        double new_i = 0.0;
        if (dist_index == 1)
        {
            new_i = (normal(rng) * post_stddev) + post_mean;
            update_residual_and_gebv(y_adj, u, col, old_i, new_i);

            chi_squared.compute(new_i * new_i);
            marker_variance(i) = chi_squared(rng);
        }
        else if (old_i != 0.0)
        {
            update_residual_and_gebv(y_adj, u, col, old_i, 0.0);
        }
        coeffs(i) = new_i;
    }

    marker_assignment.count(1) = count_1;
    marker_assignment.count(0) = static_cast<int>(coeffs.size()) - count_1;

    state.variance = detail::var(state.u)(0);
}

}  // namespace gelex::detail::Gibbs

#endif  // GELEX_ALGO_INFER_MCMC_SAMPLERS_GIBBS_B_H_
