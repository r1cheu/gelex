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

#ifndef GELEX_MODEL_BAYES_SAMPLERS_GIBBS_R_H_
#define GELEX_MODEL_BAYES_SAMPLERS_GIBBS_R_H_

#include <random>

#include "../src/model/bayes/samplers/common_op.h"
#include "../src/model/bayes/samplers/gibbs/gibbs_concept.h"
#include "../src/types/bayes_effects.h"

namespace gelex::detail::Gibbs
{

template <typename EffectT, typename StateT>
    requires IsValidEffectStatePair<EffectT, StateT>
auto R(
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
    const Eigen::VectorXd marker_variances
        = state.marker_variance(0) * effect.scale->array();
    const Eigen::Index num_components = marker_variances.size();
    Eigen::VectorXi& tracker = state.tracker;

    const auto& design_matrix = bayes::get_matrix_ref(effect.design_matrix);
    const auto& cols_norm = effect.cols_norm;

    std::normal_distribution<double> normal{0, 1};

    Eigen::VectorXd log_likelihoods(num_components);
    Eigen::VectorXd probs(num_components);
    std::vector<LikelihoodParams> likelihood_params(num_components);

    double sum_square_coeffs{};
    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        if (effect.is_monomorphic(i))
        {
            continue;
        }

        const double old_i = coeffs(i);
        const auto& col = design_matrix.col(i);

        double rhs = mkl_ddot(col, y_adj);
        if (old_i != 0.0)
        {
            rhs += cols_norm(i) * old_i;
        }

        likelihood_params[0] = {logpi(0), 0.0, 0.0};
        log_likelihoods(0) = likelihood_params[0].log_likelihood;
        for (int k = 1; k < num_components; ++k)
        {
            likelihood_params[k] = compute_likelihood_params(
                rhs,
                marker_variances(k),
                cols_norm(i),
                residual_variance,
                logpi(k));
            log_likelihoods(k) = likelihood_params[k].log_likelihood;
        }

        const double max_log_like = log_likelihoods.maxCoeff();
        probs = (log_likelihoods.array() - max_log_like).exp();

        std::discrete_distribution<int> dist(
            probs.data(), probs.data() + probs.size());
        const int dist_index = dist(rng);

        tracker(i) = dist_index;

        double new_i = 0.0;
        if (dist_index > 0)
        {
            const auto& params = likelihood_params[dist_index];

            const double post_mean = rhs * params.precision_kernel;
            const double post_stddev
                = std::sqrt(residual_variance * params.precision_kernel);

            new_i = (normal(rng) * post_stddev) + post_mean;
            update_residual_and_gebv(y_adj, u, col, old_i, new_i);
            sum_square_coeffs += (new_i * new_i) / (*effect.scale)(dist_index);
        }
        else if (old_i != 0.0)
        {
            update_residual_and_gebv(y_adj, u, col, old_i, 0.0);
        }
        coeffs(i) = new_i;
    }

    for (int k = 0; k < num_components; ++k)
    {
        state.pi.count(k) = static_cast<int>((tracker.array() == k).count());
    }

    const Eigen::Index num_nonzero = coeffs.size() - state.pi.count(0);
    detail::ScaledInvChiSq chi_squared{effect.marker_variance_prior};
    chi_squared.compute(sum_square_coeffs, num_nonzero);
    state.marker_variance(0) = chi_squared(rng);
    state.variance = detail::var(state.u)(0);
}

}  // namespace gelex::detail::Gibbs

#endif  // GELEX_MODEL_BAYES_SAMPLERS_GIBBS_R_H_
