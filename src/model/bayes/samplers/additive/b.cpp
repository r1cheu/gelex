#include "b.h"

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

auto B::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    const auto* effect = model.additive();
    auto* state = states.additive();
    auto& residual = states.residual();

    auto& y_adj = residual.y_adj;
    const double residual_variance = residual.variance;

    VectorXd logpi = state->pi.prop.array().log();

    // for convenience
    VectorXd& coeff = state->coeffs;
    auto& u = state->u;
    VectorXd& marker_variance = state->marker_variance;
    VectorXi& tracker = state->tracker;

    const auto& design_matrix = bayes::get_matrix_ref(effect->design_matrix);
    const auto& cols_norm = effect->cols_norm;

    // Setup distributions for sampling
    std::normal_distribution<double> normal{0, 1};
    std::uniform_real_distribution<double> uniform{0, 1};
    detail::ScaledInvChiSq chi_squared{effect->prior};

    for (Index i = 0; i < coeff.size(); ++i)
    {
        if (effect->is_monomorphic(i))
        {
            continue;
        }

        // for convenience
        const double old_i = coeff(i);
        const auto& col = design_matrix.col(i);
        const double variance_i = marker_variance(i);
        const double col_norm = cols_norm(i);

        const double percision_kernel
            = 1 / (col_norm + residual_variance / variance_i);

        const double post_stddev = sqrt(residual_variance * percision_kernel);
        const double logdetV
            = log((variance_i * col_norm / residual_variance) + 1);

        // calculate the posterior mean and standard deviation
        double rhs = col.dot(y_adj);
        if (old_i != 0)
        {
            rhs += col_norm * old_i;
        }
        const double post_mean = rhs * percision_kernel;

        const double L_diff
            = (-0.5 * (logdetV - post_mean * rhs / residual_variance))
              + logpi(1) - logpi(0);
        const double accept_prob = 1 / (1 + std::exp(L_diff));
        const int dist_index = (uniform(rng) < accept_prob) ? 0 : 1;
        tracker(i) = dist_index;

        double new_i = 0.0;
        if (dist_index == 1)
        {
            new_i = (normal(rng) * post_stddev) + post_mean;
            update_residual_and_gebv(y_adj, u, col, old_i, new_i);
        }
        else if (old_i != 0.0)
        {
            update_residual_and_gebv(y_adj, u, col, old_i);
        }
        coeff(i) = new_i;

        if (dist_index == 1)
        {
            chi_squared.compute(new_i * new_i);
            marker_variance(i) = chi_squared(rng);
        }
    }
    state->pi.count(1) = tracker.sum();
    state->pi.count(0) = static_cast<int>(coeff.size() - state->pi.count(1));

    state->variance = detail::var(state->u)(0);
}

}  // namespace gelex::detail::AdditiveSampler
