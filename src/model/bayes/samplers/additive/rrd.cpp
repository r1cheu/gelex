#include "rrd.h"

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

auto RRD::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    const auto* add_effect = model.additive();
    auto* add_state = states.additive();
    auto& residual = states.residual();
    const auto* dom_effect = model.dominant();
    auto* dom_state = states.dominant();

    auto& y_adj = residual.y_adj;
    const double residual_variance = residual.variance;

    VectorXd& coeff = add_state->coeffs;
    const VectorXd& dom_coeff = dom_state->coeffs;
    VectorXd& u = add_state->u;

    const double dom_ratio_mean = dom_effect->ratio_mean;
    const double dom_ratio_var = dom_state->ratio_variance;

    const double old_marker_variance = add_state->marker_variance(0);
    const auto& design_matrix
        = bayes::get_matrix_ref(add_effect->design_matrix);
    const auto& cols_norm = add_effect->cols_norm;

    std::normal_distribution<double> normal{0, 1};
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    for (Index i = 0; i < coeff.size(); ++i)
    {
        if (add_effect->is_monomorphic(i))
        {
            continue;
        }

        const double dom_i = dom_coeff(i);

        const double old_i = coeff(i);
        const auto& col = design_matrix.col(i);
        const double norm = cols_norm(i);
        const double v = norm + (residual_variance / old_marker_variance);

        // calculate the candidate posterior mean
        const double rhs = col.dot(y_adj) + (norm * old_i);
        const double post_mean = rhs / v;
        const double post_stddev = std::sqrt(residual_variance / v);

        const double cand_i = (normal(rng) * post_stddev) + post_mean;

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
    add_state->variance = detail::var(add_state->u)(0);

    detail::ScaledInvChiSq chi_squared{add_effect->prior};
    chi_squared.compute(
        add_state->coeffs.squaredNorm(),
        add_state->coeffs.size() - add_effect->num_mono());
    add_state->marker_variance(0) = chi_squared(rng);
}

}  // namespace gelex::detail::AdditiveSampler
