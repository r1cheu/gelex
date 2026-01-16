#include "gelex/estimator/freq/statistics.h"

#include <cmath>

#include <Eigen/Core>

#include "gelex/model/freq/model_new.h"
#include "gelex/optim/optimizer_state.h"

namespace gelex::statistics
{

auto compute_aic(const FreqModel& model, double loglike) -> double
{
    // k = number of variance components + number of fixed effects
    auto k = static_cast<double>(
        1 + model.random().size() + model.genetic().size()
        + model.fixed().X.cols());
    return -2.0 * loglike + 2.0 * k;
}

auto compute_bic(const FreqModel& model, double loglike) -> double
{
    auto k = static_cast<double>(
        1 + model.random().size() + model.genetic().size()
        + model.fixed().X.cols());
    auto n = static_cast<double>(model.num_individuals());
    return -2.0 * loglike + k * std::log(n);
}

auto compute_variance_se(FreqState& state, const OptimizerState& opt_state)
    -> void
{
    // se(σ) = sqrt(diag(-H⁻¹))
    // variance component order: residual, random[0..], genetic[0..]
    Eigen::VectorXd se = (-opt_state.hess_inv.diagonal()).array().sqrt();

    Eigen::Index idx = 0;

    // residual
    state.residual().variance_se = se(idx++);

    // random effects
    for (auto& r : state.random())
    {
        r.variance_se = se(idx++);
    }

    // genetic effects
    for (auto& g : state.genetic())
    {
        g.variance_se = se(idx++);
    }
}

auto compute_heritability(FreqState& state, const OptimizerState& opt_state)
    -> void
{
    // total phenotypic variance
    double sum_var = state.residual().variance;
    for (const auto& r : state.random())
    {
        sum_var += r.variance;
    }
    for (const auto& g : state.genetic())
    {
        sum_var += g.variance;
    }

    if (sum_var <= 0.0)
    {
        return;
    }

    double sum_var_sq = sum_var * sum_var;
    auto n_comp = opt_state.hess_inv.rows();

    // for each genetic effect, compute heritability and its SE using delta
    // method
    Eigen::Index genetic_start_idx
        = 1 + static_cast<Eigen::Index>(state.random().size());

    for (size_t gi = 0; gi < state.genetic().size(); ++gi)
    {
        auto& g = state.genetic()[gi];
        Eigen::Index g_idx = genetic_start_idx + static_cast<Eigen::Index>(gi);

        // h² = σ_g / Σσ
        g.heritability = g.variance / sum_var;

        // gradient for delta method:
        // ∂h²/∂σ_i = -σ_g / (Σσ)²           for i ≠ g
        // ∂h²/∂σ_g = (Σσ - σ_g) / (Σσ)²     for i == g
        Eigen::VectorXd grad(n_comp);
        for (Eigen::Index i = 0; i < n_comp; ++i)
        {
            if (i == g_idx)
            {
                grad(i) = (sum_var - g.variance) / sum_var_sq;
            }
            else
            {
                grad(i) = -g.variance / sum_var_sq;
            }
        }

        // se(h²) = sqrt(g' * (-H⁻¹) * g)
        double var_h2 = grad.dot(-opt_state.hess_inv * grad);
        g.heritability_se = std::sqrt(std::max(0.0, var_h2));
    }
}

}  // namespace gelex::statistics
