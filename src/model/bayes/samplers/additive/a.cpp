#include "a.h"

#include <random>

#include <Eigen/Core>

#include "../../bayes_effects.h"
#include "data/math_utils.h"
#include "gelex/model/bayes/model.h"
#include "model/bayes/samplers/additive/common_op.h"

namespace gelex::detail::AdditiveSampler
{
using Eigen::Index;
using Eigen::VectorXd;
using Eigen::VectorXi;

auto A::operator()(
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
    auto& u = state->u;
    VectorXd& sigma = state->marker_variance;
    const auto& design_matrix = bayes::get_matrix_ref(effect->design_matrix);
    const auto& cols_norm = effect->cols_norm;

    detail::ScaledInvChiSq chi_squared{effect->prior};
    std::normal_distribution<double> normal{0, 1};

    for (Index i = 0; i < coeff.size(); ++i)
    {
        if (effect->is_monomorphic(i))
        {
            continue;
        }
        const double old_i = coeff(i);
        const auto& col = design_matrix.col(i);
        const double col_norm = cols_norm(i);

        const double percision_kernel
            = 1 / (col_norm + residual_variance / sigma(i));

        // calculate the posterior mean and standard deviation
        const double rhs = col.dot(y_adj) + (col_norm * old_i);
        const double post_mean = rhs * percision_kernel;
        const double post_stddev = sqrt(residual_variance * percision_kernel);

        // sample a new coefficient
        const double new_i = (normal(rng) * post_stddev) + post_mean;
        coeff(i) = new_i;

        chi_squared.compute(new_i * new_i);
        sigma(i) = chi_squared(rng);
        update_residual_and_gebv(y_adj, u, col, old_i, new_i);
    }
    state->variance = detail::var(state->u)(0);
};

}  // namespace gelex::detail::AdditiveSampler
