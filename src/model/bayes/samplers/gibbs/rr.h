#ifndef GELEX_MODEL_BAYES_SAMPLERS_GIBBS_RR_H_
#define GELEX_MODEL_BAYES_SAMPLERS_GIBBS_RR_H_

#include <random>

#include "../src/model/bayes/samplers/common_op.h"
#include "../src/model/bayes/samplers/gibbs/gibbs_concept.h"
#include "../src/types/bayes_effects.h"

namespace gelex::detail::Gibbs
{

template <typename EffectT, typename StateT>
    requires IsValidEffectStatePair<EffectT, StateT>
auto RR(
    const EffectT& effect,
    StateT& state,
    bayes::ResidualState& residual,
    std::mt19937_64& rng) -> void
{
    auto& y_adj = residual.y_adj;
    const double residual_variance = residual.variance;

    Eigen::VectorXd& coeff = state.coeffs;
    const double old_marker_variance = state.marker_variance(0);
    Eigen::VectorXd& u = state.u;
    const auto& design_matrix = bayes::get_matrix_ref(effect.design_matrix);
    const auto& cols_norm = effect.cols_norm;

    const double residual_over_var = residual_variance / old_marker_variance;
    const double sqrt_residual_variance = std::sqrt(residual_variance);

    std::normal_distribution<double> normal{0, 1};

    for (Eigen::Index i = 0; i < coeff.size(); ++i)
    {
        if (effect.is_monomorphic(i))
        {
            continue;
        }

        const double old_i = coeff(i);
        const auto col = design_matrix.col(i);
        const double v = cols_norm(i) + residual_over_var;
        const double inv_v = 1.0 / v;

        const double rhs = mkl_ddot(col, y_adj) + (cols_norm(i) * old_i);
        const double post_mean = rhs * inv_v;
        const double post_stddev = sqrt_residual_variance * std::sqrt(inv_v);

        const double new_i = (normal(rng) * post_stddev) + post_mean;
        coeff(i) = new_i;
        update_residual_and_gebv(y_adj, u, col, old_i, new_i);
    }
    state.variance = detail::var(state.u)(0);

    detail::ScaledInvChiSq chi_squared{effect.marker_variance_prior};
    chi_squared.compute(coeff.squaredNorm(), coeff.size() - effect.num_mono());
    state.marker_variance(0) = chi_squared(rng);
}

}  // namespace gelex::detail::Gibbs

#endif  // GELEX_MODEL_BAYES_SAMPLERS_GIBBS_RR_H_
