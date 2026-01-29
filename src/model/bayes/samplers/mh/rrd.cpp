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

#include "../src/model/bayes/samplers/mh/rrd.h"

#include <random>

#include <Eigen/Core>

#include "../src/model/bayes/samplers/common_op.h"
#include "../src/types/bayes_effects.h"

namespace gelex::detail::MH
{
using Eigen::Index;
using Eigen::VectorXd;
using Eigen::VectorXi;

auto RRD(
    const bayes::AdditiveEffect& add_effect,
    bayes::AdditiveState& add_state,
    const bayes::DominantEffect& dom_effect,
    bayes::DominantState& dom_state,
    bayes::ResidualState& residual,
    std::mt19937_64& rng) -> void
{
    auto& y_adj = residual.y_adj;
    const double residual_variance = residual.variance;

    VectorXd& coeff = add_state.coeffs;
    const VectorXd& dom_coeff = dom_state.coeffs;
    VectorXd& u = add_state.u;
    const VectorXd& w = dom_effect.w;

    const double dom_ratio_mean = dom_effect.ratio_mean;
    const double dom_ratio_var = dom_state.ratio_variance;

    const double old_marker_variance = add_state.marker_variance(0);
    const auto& design_matrix = bayes::get_matrix_ref(add_effect.design_matrix);
    const auto col_norm = static_cast<double>(design_matrix.rows() - 1);

    std::normal_distribution<double> normal{0, 1};
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    for (Index i = 0; i < coeff.size(); ++i)
    {
        if (add_effect.is_monomorphic(i))
        {
            continue;
        }

        const double dom_i = dom_coeff(i);

        const double old_i = coeff(i);
        const auto& col = design_matrix.col(i);
        const double v = col_norm + (residual_variance / old_marker_variance);

        // calculate the candidate posterior mean
        const double rhs = blas_ddot(col, y_adj) + (col_norm * old_i);
        const double post_mean = rhs / v;
        const double post_stddev = std::sqrt(residual_variance / v);

        auto [cdf_0, pos_prob] = get_pos(w(i), dom_i, post_mean, post_stddev);
        std::bernoulli_distribution bernoulli_dist(pos_prob);

        double cand_i = 0;

        if (bernoulli_dist(rng))
        {
            std::uniform_real_distribution<double>(cdf_0, 1);
            double q = uniform(rng);
            cand_i = inverse_of_normal_cdf(q, post_mean, post_stddev);
        }
        else
        {
            std::uniform_real_distribution<double>(0, cdf_0);
            double q = uniform(rng);
            cand_i = inverse_of_normal_cdf(q, post_mean, post_stddev);
        }

        auto h = [&](double coeff_i)
        {
            if (std::abs(coeff_i) < 1e-12)
            {
                return -std::numeric_limits<double>::infinity();
            }

            const double abs_coeff = std::abs(coeff_i);
            const double coeff_sq = coeff_i * coeff_i;

            const double denom = 2 * dom_ratio_var * coeff_sq;
            const double mean_diff = dom_i - (dom_ratio_mean * abs_coeff);

            return -log(abs_coeff) - (mean_diff * mean_diff / denom);
        };

        double acceptance_ratio = std::min(1.0, std::exp(h(cand_i) - h(old_i)));
        if (uniform(rng) < acceptance_ratio)
        {
            coeff(i) = cand_i;
            update_residual_and_gebv(y_adj, u, col, old_i, cand_i);
        }
    }
    add_state.variance = detail::var(add_state.u)(0);

    detail::ScaledInvChiSq chi_squared{add_effect.marker_variance_prior};
    chi_squared.compute(
        add_state.coeffs.squaredNorm(),
        add_state.coeffs.size() - add_effect.num_mono());
    add_state.marker_variance(0) = chi_squared(rng);
}

}  // namespace gelex::detail::MH
