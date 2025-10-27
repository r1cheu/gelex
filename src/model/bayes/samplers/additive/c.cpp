#include "c.h"

#include <Eigen/Core>
#include <cmath>
#include <random>
#include <ranges>

#include "../../bayes_effects.h"
#include "../src/utils/math_utils.h"
#include "gelex/model/bayes/model.h"
#include "model/bayes/samplers/additive/common_op.h"

namespace gelex::detail::AdditiveSampler
{
using Eigen::Index;
using Eigen::VectorXd;
using Eigen::VectorXi;

auto C::operator()(
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
    const double marker_variance = state->marker_variance(0);
    VectorXi& tracker = state->tracker;

    const auto& design_matrix = bayes::get_matrix_ref(effect->design_matrix);
    const auto col_norm = static_cast<double>(design_matrix.rows() - 1);

    std::normal_distribution<double> normal{0, 1};
    std::uniform_real_distribution<double> uniform{0, 1};

    const double residual_over_marker_variance
        = residual_variance / marker_variance;

    double sum_square_coeffs{};

    for (const Index i : std::views::iota(0, coeffs.size()))
    {
        if (effect->is_monomorphic(i))
        {
            continue;
        }

        const double old_i = coeffs(i);
        const auto& col = design_matrix.col(i);

        double rhs = col.dot(y_adj);
        if (old_i != 0.0)
        {
            rhs += col_norm * old_i;
        }

        auto [post_mean, post_stddev, log_like_kernel]
            = compute_posterior_params_core(
                rhs,
                col_norm,
                residual_variance,
                residual_over_marker_variance);

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
            sum_square_coeffs += new_i * new_i;
        }
        else if (old_i != 0.0)
        {
            update_residual_and_gebv(y_adj, u, col, old_i, 0.0);
        }
        coeffs(i) = new_i;
    }

    state->pi.count(1) = tracker.sum();
    state->pi.count(0) = static_cast<int>(coeffs.size() - state->pi.count(1));

    detail::ScaledInvChiSq chi_squared{effect->prior};
    chi_squared.compute(sum_square_coeffs, state->pi.count(1));
    state->marker_variance(0) = chi_squared(rng);

    state->variance = detail::var(state->u)(0);
}

}  // namespace gelex::detail::AdditiveSampler
