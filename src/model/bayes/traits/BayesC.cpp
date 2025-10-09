#include "BayesC.h"
#include <fmt/format.h>
#include <Eigen/Core>

#include "../src/logger/logger_utils.h"
#include "../src/model/bayes/bayes_effects.h"
#include "../src/model/bayes/distribution.h"
#include "data/math_utils.h"

namespace gelex
{
using Eigen::Ref;

using Eigen::Index;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

VectorXd BayesCTrait::default_marker_variance(Index n_snp) const
{
    VectorXd sigma = VectorXd::Zero(1);
    return sigma;
}

VectorXd BayesCTrait::default_pi() const
{
    VectorXd pi{{0.95, 0.05}};
    return pi;
}

bool BayesCTrait::estimate_pi() const
{
    return false;
}

std::vector<std::string> BayesCTrait::prior_info(
    double nu,
    double s2,
    const Ref<const VectorXd>& pi) const
{
    std::vector<std::string> prior_strings;
    prior_strings.reserve(3);
    prior_strings.emplace_back("BayesC");
    prior_strings.emplace_back("      ├─ αᵢ ~ 0.05 N(0, σ²) + 0.95 δ₀");
    prior_strings.emplace_back(
        fmt::format("      └─ {}", detail::sigma_prior("", nu, s2)));
    return prior_strings;
}

void BayesCTrait::operator()(
    const bayes::AdditiveEffect& effect,
    bayes::AdditiveStatus& state,
    Ref<VectorXd> y_adj,
    double sigma_e,
    std::mt19937_64& rng) const
{
    VectorXd logpi = state.pi.prop.array().log();

    // for convenience
    VectorXd& coeff = state.coeff;
    auto& u = state.u;
    const double sigma_g = state.marker_variance(0);
    VectorXi& tracker = state.tracker;

    const auto& design_matrix = effect.design_matrix.matrix();

    // col_norm is n - 1, since we have normalized the design matrix
    const auto& cols_norm = effect.cols_norm;
    const VectorXd percision_kernel
        = 1.0 / (cols_norm.array() + sigma_e / sigma_g);

    // calculate the posterior standard deviation
    const VectorXd post_stddev = (sigma_e * percision_kernel.array()).sqrt();
    const VectorXd logdetV
        = ((sigma_g * cols_norm.array() / sigma_e) + 1).log();

    // Setup distributions for sampling
    std::normal_distribution<double> normal{0, 1};
    std::uniform_real_distribution<double> uniform{0, 1};

    double var_g{};
    for (Index i = 0; i < coeff.size(); ++i)
    {
        if (effect.is_monomorphic(i))
        {
            continue;
        }

        // for convenience
        const double old_i = coeff(i);
        const auto& col = design_matrix.col(i);
        const double col_norm = cols_norm(i);

        // calculate the posterior mean
        double rhs = col.dot(y_adj);
        if (old_i != 0)
        {
            rhs += col_norm * old_i;
        }
        const double post_mean = rhs * percision_kernel(i);

        // determine which distribution to sample from
        const double L_diff = (-0.5 * (logdetV(i) - post_mean * rhs / sigma_e))
                              + logpi(1) - logpi(0);
        const double accept_prob = 1 / (1 + std::exp(L_diff));
        const int dist_index = (uniform(rng) < accept_prob) ? 0 : 1;
        tracker(i) = dist_index;

        // sample a new coefficient
        double new_i = 0.0;
        if (dist_index == 1)
        {
            new_i = (normal(rng) * post_stddev(i)) + post_mean;
            const double diff = old_i - new_i;
            y_adj.array() += col.array() * diff;
            u.array() -= col.array() * diff;
            var_g += new_i * new_i;
        }
        else if (old_i != 0.0)
        {
            y_adj.array() += col.array() * old_i;
            u.array() -= col.array() * old_i;
        }
        coeff(i) = new_i;
    }
    state.pi.count(1) = tracker.sum();
    state.pi.count(0) = static_cast<int>(coeff.size() - state.pi.count(1));

    // sample variance
    detail::ScaledInvChiSq chi_squared{effect.prior};
    chi_squared.compute(var_g, state.pi.count(1));
    state.marker_variance(0) = chi_squared(rng);

    state.effect_variance = detail::var(state.u)(0);
}

}  // namespace gelex
