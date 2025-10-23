#include "common.h"

#include <random>

#include <Eigen/Core>

#include "../bayes_effects.h"
#include "gelex/model/bayes/model.h"

namespace gelex::detail::CommonSampler
{

using Eigen::Index;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

auto Fixed::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    if (const auto* effect = model.fixed(); effect == nullptr)
    {
        return;
    }

    const auto* effect = model.fixed();
    auto* state = states.fixed();
    auto& residual = states.residual();

    auto& y_adj = residual.y_adj;
    const double residual_variance = residual.variance;

    auto& coeffs = state->coeffs;
    const auto& cols_norm = effect->cols_norm;
    const auto& design_matrix = effect->design_matrix;

    std::normal_distribution<double> normal{0, 1};

    for (int i = 0; i < coeffs.size(); ++i)
    {
        const double old_i = coeffs(i);
        const auto& col = design_matrix.col(i);
        const double norm = cols_norm(i);

        const double rhs = col.dot(y_adj) + (norm * old_i);
        const double post_mean = rhs / norm;
        const double post_stddev = std::sqrt(residual_variance / norm);

        // sample a new coefficient
        const double new_i = (normal(rng) * post_stddev) + post_mean;
        coeffs(i) = new_i;

        // update the y_adj vector
        const double diff = old_i - new_i;
        y_adj.array() += diff * col.array();
    }
}

auto Random::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    if (const auto& effect = model.random(); effect.empty())
    {
        return;
    }

    const auto& effects = model.random();
    auto& state = states.random();
    auto& residual = states.residual();

    for (Index i = 0; i < static_cast<Index>(effects.size()); ++i)
    {
        sample_impl(effects[i], state[i], residual, rng);
    }
}

auto Random::sample_impl(
    const bayes::RandomEffect& effect,
    bayes::RandomState& status,
    bayes::ResidualState& residual,
    std::mt19937_64& rng) -> void
{
    VectorXd& coeff = status.coeffs;
    const VectorXd& cols_norm = effect.cols_norm;
    const MatrixXd& design_matrix = effect.design_matrix;

    VectorXd& y_adj = residual.y_adj;
    const double residual_variance = residual.variance;

    const double sigma = status.variance;

    // calculate precision kernel and posterior standard deviation
    const VectorXd inv_scaler
        = 1.0 / (cols_norm.array() + residual_variance / sigma);
    const VectorXd post_stddev
        = (residual_variance * inv_scaler.array()).sqrt();

    // Setup distributions for sampling
    std::normal_distribution<double> normal{0, 1};

    for (int i = 0; i < coeff.size(); ++i)
    {
        // for convenience
        const double old_i = coeff(i);
        const auto& col = design_matrix.col(i);
        const double norm = cols_norm(i);

        // calculate the posterior mean
        const double rhs = col.dot(y_adj) + (norm * old_i);
        const double post_mean = rhs * inv_scaler(i);

        // sample a new coefficient
        const double new_i = (normal(rng) * post_stddev(i)) + post_mean;
        coeff(i) = new_i;

        // update the y_adj vector
        const double diff = old_i - new_i;
        y_adj.array() += col.array() * diff;
    }

    detail::ScaledInvChiSq chi_squared{effect.prior};
    chi_squared.compute(coeff.squaredNorm(), coeff.size());
    status.variance = chi_squared(rng);
}

auto Residual::operator()(
    const BayesModel& model,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    auto& residual = states.residual();
    detail::ScaledInvChiSq chi_squared{model.residual().prior};
    chi_squared.compute(residual.y_adj.squaredNorm(), model.num_individuals());
    residual.variance = chi_squared(rng);
}
}  // namespace gelex::detail::CommonSampler
