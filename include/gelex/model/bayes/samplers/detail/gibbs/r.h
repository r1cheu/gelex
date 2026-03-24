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

#ifndef GELEX_MODEL_BAYES_SAMPLERS_DETAIL_GIBBS_R_H_
#define GELEX_MODEL_BAYES_SAMPLERS_DETAIL_GIBBS_R_H_

#include <array>
#include <random>
#include <type_traits>

#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/samplers/detail/common_op.h"
#include "gelex/model/bayes/samplers/detail/gibbs/gibbs_concept.h"
#include "gelex/model/bayes/samplers/detail/gibbs/r_policy.h"
#include "gelex/model/bayes/states.h"

namespace gelex::detail::Gibbs
{

template <
    typename Policy = policy::Symmetric,
    typename EffectT,
    typename StateT>
    requires IsValidEffectStatePair<EffectT, StateT>
auto R(
    const EffectT& effect,
    const bayes::GeneticPrior& prior,
    StateT& state,
    bayes::ResidualState& residual,
    std::mt19937_64& rng) -> void
{
    auto& y_adj = residual.y_adj;
    const double residual_variance = residual.variance;

    auto& marker_assignment
        = std::get<bayes::ComponentAllocation>(*state.group);
    Eigen::VectorXd& coeffs = state.coeffs;
    auto& u = state.u;
    const auto& multiplier
        = std::get<bayes::MixturePrior>(prior.marker).multiplier;
    const Eigen::VectorXd marker_variances
        = state.marker_variance(0) * multiplier.array();
    const Eigen::Index num_scale_classes = marker_variances.size();
    auto& tracker = marker_assignment.assignment.tracker;

    const auto& X = bayes::get_matrix_ref(effect.X);
    const auto& cols_squared_norm = effect.cols_squared_norm;

    std::uniform_real_distribution<double> uniform{0.0, 1.0};

    Policy policy;
    if constexpr (std::is_same_v<Policy, policy::AT>)
    {
        policy.init(*state.sign);
    }
    else
    {
        policy.init();
    }
    const int num_options
        = Policy::num_options(static_cast<int>(num_scale_classes));

    std::array<PosteriorParams, kMaxMixtureComponents> scale_pp{};
    std::array<double, kMaxMixtureComponents> logpi_buf{};
    std::array<double, kMaxMixtureComponents> ll_buf{};
    std::array<double, kMaxMixtureComponents> probs_buf{};
    std::array<int, kMaxMixtureComponents> pi_count{};

    for (Eigen::Index k = 0; k < num_scale_classes; ++k)
    {
        logpi_buf[k] = std::log(marker_assignment.assignment.proportion(k));
    }

    double sum_square_coeffs{};
    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        if (effect.is_monomorphic(i))
        {
            pi_count[tracker(i)]++;
            continue;
        }

        const double old_i = coeffs(i);
        const auto& col = X.col(i);

        double rhs = blas_ddot(col, y_adj);
        if (old_i != 0.0)
        {
            rhs += cols_squared_norm(i) * old_i;
        }

        for (int c = 1; c < num_scale_classes; ++c)
        {
            scale_pp[c] = compute_posterior_params(
                rhs,
                marker_variances(c),
                cols_squared_norm(i),
                residual_variance);
        }

        double max_ll = policy.compute_log_probs(
            ll_buf, logpi_buf, scale_pp, static_cast<int>(num_scale_classes));

        double total = 0.0;
        for (int k = 0; k < num_options; ++k)
        {
            probs_buf[k] = std::exp(ll_buf[k] - max_ll);
            total += probs_buf[k];
        }

        const double u_val = uniform(rng) * total;
        int8_t dist_index = num_options - 1;
        double cumsum = 0.0;
        for (int k = 0; k < num_options; ++k)
        {
            cumsum += probs_buf[k];
            if (u_val < cumsum)
            {
                dist_index = k;
                break;
            }
        }

        const int8_t old_index = tracker(i);

        auto [class_index, new_i]
            = (dist_index > 0) ? policy.sample_effect(dist_index, scale_pp, rng)
                               : std::pair<int, double>{0, 0.0};

        tracker(i) = class_index;
        pi_count[class_index]++;

        if (dist_index > 0)
        {
            update_residual_and_gebv(y_adj, u, col, old_i, new_i);
            sum_square_coeffs += (new_i * new_i) / multiplier(class_index);
        }
        else if (old_i != 0.0)
        {
            update_residual_and_gebv(y_adj, u, col, old_i, 0.0);
        }
        coeffs(i) = new_i;

        if constexpr (std::is_same_v<Policy, policy::AT>)
        {
            int8_t new_sign_index{};
            if (dist_index == 0)
            {
                new_sign_index = 0;
            }
            else if (dist_index % 2 == 1)
            {
                new_sign_index = 1;
            }
            else
            {
                new_sign_index = 2;
            }
            state.sign->tracker(i) = new_sign_index;
        }

        update_component_u(
            marker_assignment.component_u,
            old_index,
            old_i,
            class_index,
            new_i,
            col);
    }

    for (int k = 0; k < num_scale_classes; ++k)
    {
        marker_assignment.assignment.count(k) = pi_count[k];
    }

    const Eigen::Index num_nonzero
        = coeffs.size() - marker_assignment.assignment.count(0);
    const auto& marker_prior = std::get<bayes::MixturePrior>(prior.marker);
    detail::ScaledInvChiSq chi_squared{marker_prior.variance.param};
    chi_squared.compute(sum_square_coeffs, num_nonzero);
    state.marker_variance(0) = chi_squared(rng);
    state.variance = detail::var(state.u)(0);

    compute_component_variances(state);

    if constexpr (std::is_same_v<Policy, policy::AT>)
    {
        policy.finalize(*state.sign, rng);
    }
    else
    {
        policy.finalize(rng);
    }
}

}  // namespace gelex::detail::Gibbs

#endif  // GELEX_MODEL_BAYES_SAMPLERS_DETAIL_GIBBS_R_H_
