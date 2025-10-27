#include "b.h"

#include <Eigen/Core>
#include <cmath>  // std::log, std::sqrt, std::exp
#include <random>
#include <ranges>  // std::views::iota

#include "../../bayes_effects.h"
#include "../src/utils/math_utils.h"
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

    const VectorXd logpi = state->pi.prop.array().log();

    VectorXd& coeffs = state->coeffs;
    auto& u = state->u;
    VectorXd& marker_variance = state->marker_variance;
    VectorXi& tracker = state->tracker;

    const auto& design_matrix = bayes::get_matrix_ref(effect->design_matrix);
    const auto& cols_norm = effect->cols_norm;

    std::normal_distribution<double> normal{0, 1};
    std::uniform_real_distribution<double> uniform{0, 1};
    detail::ScaledInvChiSq chi_squared{effect->prior};

    for (const Index i : std::views::iota(0, coeffs.size()))
    {
        if (effect->is_monomorphic(i))
        {
            continue;
        }

        const double old_i = coeffs(i);
        const auto& col = design_matrix.col(i);
        const double variance_i = marker_variance(i);
        const double col_norm = cols_norm(i);

        double rhs = col.dot(y_adj);
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

    state->pi.count(1) = tracker.sum();
    state->pi.count(0) = static_cast<int>(coeffs.size() - state->pi.count(1));

    state->variance = detail::var(state->u)(0);
}

}  // namespace gelex::detail::AdditiveSampler
