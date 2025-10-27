#include "r.h"

#include <random>
#include <ranges>

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

auto R::operator()(
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
    const VectorXd marker_variances
        = state->marker_variance(0) * effect->scale.array();
    const Index num_components = marker_variances.size();
    VectorXi& tracker = state->tracker;

    const auto& design_matrix = bayes::get_matrix_ref(effect->design_matrix);
    const auto& cols_norm = effect->cols_norm;

    std::normal_distribution<double> normal{0, 1};

    VectorXd log_likelihoods(num_components);
    VectorXd probs(num_components);
    std::vector<LikelihoodParams> likelihood_params(num_components);

    double sum_square_coeffs{};
    for (const Index i : std::views::iota(0, coeffs.size()))
    {
        if (effect->is_monomorphic(i))
        {
            continue;
        }

        const double old_i = coeffs(i);
        const auto& col = design_matrix.col(i);
        const double col_norm = cols_norm(i);

        double rhs = col.dot(y_adj);
        if (old_i != 0.0)
        {
            rhs += col_norm * old_i;
        }

        likelihood_params[0] = {logpi(0), 0.0, 0.0};
        log_likelihoods(0) = likelihood_params[0].log_likelihood;
        for (int k = 1; k < num_components; ++k)
        {
            likelihood_params[k] = compute_likelihood_params(
                rhs,
                marker_variances(k),
                col_norm,
                residual_variance,
                logpi(k));
            log_likelihoods(k) = likelihood_params[k].log_likelihood;
        }

        const double max_log_like = log_likelihoods.maxCoeff();
        probs = (log_likelihoods.array() - max_log_like).exp();

        std::discrete_distribution<int> dist(
            probs.data(), probs.data() + probs.size());
        const int dist_index = dist(rng);

        tracker(i) = dist_index;

        double new_i = 0.0;
        if (dist_index > 0)
        {
            const auto& params = likelihood_params[dist_index];

            const double post_mean = rhs * params.precision_kernel;
            const double post_stddev
                = std::sqrt(residual_variance * params.precision_kernel);

            new_i = (normal(rng) * post_stddev) + post_mean;
            update_residual_and_gebv(y_adj, u, col, old_i, new_i);
            sum_square_coeffs += (new_i * new_i) / effect->scale(dist_index);
        }
        else if (old_i != 0.0)
        {
            update_residual_and_gebv(y_adj, u, col, old_i, 0.0);
        }
        coeffs(i) = new_i;
    }

    for (int k = 0; k < num_components; ++k)
    {
        state->pi.count(k) = static_cast<int>((tracker.array() == k).count());
    }

    const Index num_nonzero = coeffs.size() - state->pi.count(0);
    detail::ScaledInvChiSq chi_squared{effect->prior};
    chi_squared.compute(sum_square_coeffs, num_nonzero);
    state->marker_variance(0) = chi_squared(rng);

    state->variance = detail::var(state->u)(0);
}

}  // namespace gelex::detail::AdditiveSampler
