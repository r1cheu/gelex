#pragma once

#include <random>

#include "gibbs_concept.h"
#include "model/bayes/samplers/common_op.h"
#include "types/bayes_effects.h"

namespace gelex::detail::Gibbs
{

template <typename EffectT, typename StateT>
    requires IsValidEffectStatePair<EffectT, StateT>
auto B(
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
    Eigen::VectorXd& marker_variance = state.marker_variance;
    Eigen::VectorXi& tracker = state.tracker;

    const auto& design_matrix = bayes::get_matrix_ref(effect.design_matrix);
    const auto col_norm = static_cast<double>(design_matrix.rows() - 1);

    std::normal_distribution<double> normal{0, 1};
    std::uniform_real_distribution<double> uniform{0, 1};
    detail::ScaledInvChiSq chi_squared{effect.marker_variance_prior};

    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        if (effect.is_monomorphic(i))
        {
            continue;
        }

        const double old_i = coeffs(i);
        const auto& col = design_matrix.col(i);
        const double variance_i = marker_variance(i);

        double rhs = mkl_ddot(col, y_adj);
        if (old_i != 0.0)
        {
            rhs += col_norm * old_i;
        }

        auto [post_mean, post_stddev, log_like_kernel]
            = compute_posterior_params(
                rhs, variance_i, col_norm, residual_variance);

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

            chi_squared.compute(new_i * new_i);
            marker_variance(i) = chi_squared(rng);
        }
        else if (old_i != 0.0)
        {
            update_residual_and_gebv(y_adj, u, col, old_i, 0.0);
        }
        coeffs(i) = new_i;
    }

    state.pi.count(1) = tracker.sum();
    state.pi.count(0) = static_cast<int>(coeffs.size() - state.pi.count(1));

    state.variance = detail::var(state.u)(0);
}

}  // namespace gelex::detail::Gibbs
