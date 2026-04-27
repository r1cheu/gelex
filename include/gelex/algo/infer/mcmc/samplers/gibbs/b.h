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

#include "gelex/algo/infer/detail/marker_op.h"
#include "gelex/algo/infer/mcmc/samplers/gibbs/gibbs_concept.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/infra/stats/descriptive.h"
#include "gelex/model/bayes/genotype_storage.h"
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

    std::uniform_real_distribution<double> uniform{0, 1};
    const auto& marker_prior = std::get<bayes::SpikePrior>(prior.marker);
    NormalSampler<double> normal_sampler(0.0);
    ScaledInvChi2Sampler<double> variance_sampler(
        marker_prior.variance.param.nu, marker_prior.variance.param.s2);

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

        double rhs = blas_ddot(col, y_adj);
        if (old_i != 0.0)
        {
            rhs += cols_squared_norm(i) * old_i;
        }

        const auto post
            = normal_sampler.set_prior_var(marker_variance(i))
                  .posterior_with_logL(
                      {cols_squared_norm(i), rhs, residual_variance});

        const double log_like_1_minus_0
            = post.log_likelihood_kernel + logpi_1 - logpi_0;
        const double prob_component_0
            = 1.0 / (1.0 + std::exp(log_like_1_minus_0));

        const int8_t dist_index = (uniform(rng) < prob_component_0) ? 0 : 1;
        tracker(i) = dist_index;
        count_1 += dist_index;

        double new_i = 0.0;
        if (dist_index == 1)
        {
            new_i = normal_sampler.draw(post.params, rng);
            update_residual_and_gebv(y_adj, u, col, old_i, new_i);
            marker_variance(i) = variance_sampler({1, new_i * new_i}, rng);
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
