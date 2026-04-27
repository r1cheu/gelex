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

#ifndef GELEX_ALGO_INFER_MCMC_SAMPLERS_GIBBS_RR_H_
#define GELEX_ALGO_INFER_MCMC_SAMPLERS_GIBBS_RR_H_

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
auto RR(
    const EffectT& effect,
    const bayes::GeneticPrior& prior,
    StateT& state,
    bayes::ResidualState& residual,
    std::mt19937_64& rng) -> void
{
    auto& y_adj = residual.y_adj;
    const double residual_variance = residual.variance;

    Eigen::VectorXd& coeff = state.coeffs;
    Eigen::VectorXd& u = state.u;
    const auto& X = bayes::get_matrix_ref(effect.X);
    const auto& cols_squared_norm = effect.cols_squared_norm;

    NormalSampler<double> normal_sampler(state.marker_variance(0));

    for (Eigen::Index i = 0; i < coeff.size(); ++i)
    {
        if (effect.is_monomorphic(i))
        {
            continue;
        }

        const double old_i = coeff(i);
        const auto col = X.col(i);

        const double rhs
            = blas_ddot(col, y_adj) + (cols_squared_norm(i) * old_i);
        const double new_i = normal_sampler(
            {cols_squared_norm(i), rhs, residual_variance}, rng);
        coeff(i) = new_i;
        update_residual_and_gebv(y_adj, u, col, old_i, new_i);
    }
    state.variance = detail::var(state.u)(0);

    const auto& marker_prior = std::get<bayes::ContinuousPrior>(prior.marker);
    ScaledInvChi2Sampler<double> variance_sampler(
        marker_prior.variance.param.nu, marker_prior.variance.param.s2);
    state.marker_variance(0) = variance_sampler(
        {coeff.size() - effect.num_mono(), coeff.squaredNorm()}, rng);
}

}  // namespace gelex::detail::Gibbs

#endif  // GELEX_ALGO_INFER_MCMC_SAMPLERS_GIBBS_RR_H_
