#include "dominant.h"
#include <random>
#include <ranges>

#include <Eigen/Core>

#include "../src/utils/math_utils.h"
#include "gelex/model/bayes/model.h"
#include "model/bayes/samplers/common_op.h"
#include "types/bayes_effects.h"

#include "model/bayes/samplers/gibbs/a.h"
#include "model/bayes/samplers/gibbs/b.h"
#include "model/bayes/samplers/gibbs/c.h"
#include "model/bayes/samplers/gibbs/r.h"
#include "model/bayes/samplers/gibbs/rr.h"

namespace gelex::detail::DominantSampler
{

using Eigen::Index;
using Eigen::VectorXd;
using Eigen::VectorXi;

auto A::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    const auto* effect = model.dominant();
    auto* state = states.dominant();
    auto& residual = states.residual();

    Gibbs::A(*effect, *state, residual, rng);

    // Update ratios using optimized vector operation
}

auto B::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    const auto* effect = model.dominant();
    auto* state = states.dominant();
    auto& residual = states.residual();

    Gibbs::B(*effect, *state, residual, rng);

    // Update ratios using optimized vector operation
}

auto C::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    const auto* effect = model.dominant();
    auto* state = states.dominant();
    auto& residual = states.residual();

    Gibbs::C(*effect, *state, residual, rng);
}

auto R::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    const auto* effect = model.dominant();
    auto* state = states.dominant();
    auto& residual = states.residual();

    Gibbs::R(*effect, *state, residual, rng);
}

auto RR::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    const auto* effect = model.dominant();
    auto* state = states.dominant();
    auto& residual = states.residual();

    Gibbs::RR(*effect, *state, residual, rng);
}

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
    VectorXd& u = dom_state->u;
    const VectorXd& add_coeffs = add_state->coeffs;
    const VectorXd& w = dom_effect->w;

    const double ratio_mean = dom_effect->ratio_mean;
    const double ratio_var = dom_state->ratio_variance;

    const auto& design_matrix
        = bayes::get_matrix_ref(dom_effect->design_matrix);
    const auto col_norm = static_cast<double>(design_matrix.rows() - 1);
    auto& ratios = dom_state->ratios;

    std::normal_distribution<double> normal{0, 1};
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    for (Index i : std::views::iota(0, coeffs.size()))
    {
        if (dom_effect->is_monomorphic(i))
        {
            continue;
        }

        const double old_i = coeffs(i);
        const double add_i = add_coeffs(i);

        // don't sample dominant if the additive effect is zero
        if (std::abs(add_i) < 1e-12)
        {
            coeffs(i) = 0.0;
            const double diff = old_i - 0.0;
            if (std::abs(diff) > 1e-12)
            {
                update_residual_and_gebv(
                    y_adj, u, design_matrix.col(i), old_i, 0.0);
            }
            continue;
        }

        const auto& col = design_matrix.col(i);

        const double residual_over_var
            = residual_variance / (ratio_var * add_i * add_i);
        const double v = col_norm + residual_over_var;

        const double rhs = (col.dot(y_adj)) + (col_norm * old_i)
                           + (std::abs(add_i) * ratio_mean * residual_over_var);

        const double post_mean = rhs / v;
        const double post_stddev = std::sqrt(residual_variance / v);
        auto [cdf_0, pos_prob] = get_pos(w(i), add_i, post_mean, post_stddev);
        std::bernoulli_distribution bernoulli_dist(pos_prob);

        double new_i = 0;

        if (bernoulli_dist(rng))
        {
            std::uniform_real_distribution<double>(cdf_0, 1);
            double q = uniform(rng);
            new_i = inverse_of_normal_cdf(q, post_mean, post_stddev);
        }
        else
        {
            std::uniform_real_distribution<double>(0, cdf_0);
            double q = uniform(rng);
            new_i = inverse_of_normal_cdf(q, post_mean, post_stddev);
        }

        coeffs(i) = new_i;
        update_residual_and_gebv(y_adj, u, col, old_i, new_i);
    }

    // Update all ratios using optimized vector operation
    ratios = detail::compute_dominant_ratios(coeffs, add_coeffs);

    dom_state->variance = detail::var(dom_state->u)(0);
}

auto RatioMean::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    const auto* effect = model.dominant();
    auto* state = states.dominant();

    auto num_coeffs = static_cast<double>(
        bayes::get_cols(effect->design_matrix)
        - bayes::num_mono_variant(effect->design_matrix));
    const double ratio_var = state->ratio_variance;

    const auto& mean_prior = effect->mean_prior;

    const double post_stddev
        = std::sqrt(ratio_var / (ratio_var / mean_prior.var + num_coeffs));

    const double post_mean
        = (mean_prior.mean * (ratio_var / mean_prior.var) + state->ratios.sum())
          / (ratio_var / mean_prior.var + num_coeffs);

    std::normal_distribution<double> normal{0, 1};
    state->ratio_mean = post_mean + post_stddev * normal(rng);
}

auto RatioVar::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    const auto* effect = model.dominant();
    auto* state = states.dominant();

    detail::ScaledInvChiSq dist(effect->var_prior);

    const double sum_of_squared_errors
        = (state->ratios.array() - state->ratio_mean).square().sum();

    const auto num_coeffs
        = (bayes::get_cols(effect->design_matrix)
           - bayes::num_mono_variant(effect->design_matrix));
    dist.compute(sum_of_squared_errors, num_coeffs);
    state->ratio_variance = dist(rng);
}
}  // namespace gelex::detail::DominantSampler
