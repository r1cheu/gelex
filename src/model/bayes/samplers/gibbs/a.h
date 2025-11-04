#pragma once

#include <random>

#include "gibbs_concept.h"
#include "model/bayes/samplers/common_op.h"
#include "types/bayes_effects.h"

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
    const auto& design_matrix = bayes::get_matrix_ref(effect.design_matrix);
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
        const auto& col = design_matrix.col(i);

        const double percision_kernel
            = 1 / (cols_norm(i) + residual_variance / sigma(i));

        // calculate the posterior mean and standard deviation
        const double rhs = mkl_ddot(col, y_adj) + (cols_norm(i) * old_i);
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
