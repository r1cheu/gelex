#include "rr.h"

#include <random>

#include <Eigen/Core>

#include "../../bayes_effects.h"
#include "../src/utils/math_utils.h"
#include "gelex/model/bayes/model.h"
#include "model/bayes/samplers/additive/common_op.h"

namespace gelex::detail::AdditiveSampler
{
using Eigen::Index;
using Eigen::VectorXd;
using Eigen::VectorXi;

auto RR::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    const auto* effect = model.additive();
    auto* state = states.additive();
    auto& residual = states.residual();

    auto& y_adj = residual.y_adj;
    const double residual_variance = residual.variance;

    VectorXd& coeff = state->coeffs;
    const double old_marker_variance = state->marker_variance(0);
    VectorXd& u = state->u;
    const auto& design_matrix = bayes::get_matrix_ref(effect->design_matrix);
    const auto col_norm = static_cast<double>(design_matrix.rows() - 1);

    const double residual_over_var = residual_variance / old_marker_variance;
    const double sqrt_residual_variance = std::sqrt(residual_variance);

    std::normal_distribution<double> normal{0, 1};

    for (Index i = 0; i < coeff.size(); ++i)
    {
        if (effect->is_monomorphic(i))
        {
            continue;
        }

        const double old_i = coeff(i);
        const auto col = design_matrix.col(i);
        const double v = col_norm + residual_over_var;
        const double inv_v = 1.0 / v;

        const double rhs = col.dot(y_adj) + (col_norm * old_i);
        const double post_mean = rhs * inv_v;
        const double post_stddev = sqrt_residual_variance * std::sqrt(inv_v);

        const double new_i = (normal(rng) * post_stddev) + post_mean;
        coeff(i) = new_i;
        update_residual_and_gebv(y_adj, u, col, old_i, new_i);
    }
    state->variance = detail::var(state->u)(0);

    detail::ScaledInvChiSq chi_squared{effect->prior};
    chi_squared.compute(coeff.squaredNorm(), coeff.size() - effect->num_mono());
    state->marker_variance(0) = chi_squared(rng);
}

}  // namespace gelex::detail::AdditiveSampler
