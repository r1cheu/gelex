#include "dominant.h"
#include <random>
#include <ranges>

#include <Eigen/Core>

#include "../bayes_effects.h"
#include "data/math_utils.h"
#include "gelex/model/bayes/model.h"

namespace gelex::detail::DominantSampler
{

using Eigen::Index;
using Eigen::VectorXd;
using Eigen::VectorXi;

auto Coeff::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    const auto* dom_effect = model.dominant();
    auto* dom_state = states.dominant();
    const auto* add_state = states.additive();

    auto& residual = states.residual();
    auto& y_adj = residual.y_adj;
    const double residual_variance = residual.variance;

    VectorXd& coeffs = dom_state->coeffs;
    VectorXd& dom_ratio = dom_state->ratios;
    VectorXd& u = dom_state->u;
    const VectorXd& add_coeff = add_state->coeffs;

    const double ratio_percision = 1 / dom_state->ratio_variance;
    const double e_over_ratio_var = residual_variance * ratio_percision;

    const auto& design_matrix
        = bayes::get_matrix_ref(dom_effect->design_matrix);

    const VectorXd scaled_cols_norm
        = dom_effect->cols_norm.array() * add_coeff.array().square();

    std::normal_distribution<double> normal{0, 1};

    for (Index i : std::views::iota(0, coeffs.size()))
    {
        if (dom_effect->is_monomorphic(i))
        {
            continue;
        }

        const double old_i = coeffs(i);
        const double old_add_i = add_coeff(i);

        // don't sample dominant if the additive effect is zero
        if (std::abs(old_add_i) < 1e-12)
        {
            dom_ratio(i) = 0.0;
            coeffs(i) = 0.0;

            const double diff = old_i - 0.0;
            if (std::abs(diff) > 1e-12)
            {
                const auto& col = design_matrix.col(i);
                y_adj.array() += col.array() * diff;
                u.array() -= col.array() * diff;
            }
            continue;
        }

        const double old_ratio_i = dom_ratio(i);
        const auto& col = design_matrix.col(i);
        const double norm = scaled_cols_norm(i);

        const double v = norm + e_over_ratio_var;

        const double rhs = (std::abs(old_add_i) * col.dot(y_adj))
                           + (norm * old_ratio_i)
                           + (dom_state->ratio_mean * e_over_ratio_var);

        const double post_mean = rhs / v;
        const double post_stddev = std::sqrt(residual_variance / v);

        const double new_ratio_i = (normal(rng) * post_stddev) + post_mean;
        const double new_i = new_ratio_i * std::abs(old_add_i);
        dom_ratio(i) = new_ratio_i;
        coeffs(i) = new_i;

        const double diff = old_i - new_i;
        if (std::abs(diff) > 1e-12)
        {
            y_adj.array() += col.array() * diff;
            u.array() -= col.array() * diff;
        }
    }
    dom_state->variance = detail::var(dom_state->u)(0);
}
}  // namespace gelex::detail::DominantSampler
