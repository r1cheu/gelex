#include "BayesB.h"

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

VectorXd BayesBTrait::default_sigma(Index n_snp) const
{
    VectorXd sigma = VectorXd::Zero(n_snp);
    return sigma;
}

VectorXd BayesBTrait::default_pi() const
{
    VectorXd pi{{0.95, 0.05}};
    return pi;
}

bool BayesBTrait::estimate_pi() const
{
    return false;
}

std::vector<std::string> BayesBTrait::prior_info(
    double nu,
    double s2,
    const Ref<const VectorXd>& pi) const
{
    std::vector<std::string> prior_strings;
    prior_strings.reserve(3);
    prior_strings.emplace_back("BayesB");
    prior_strings.emplace_back("      ├─ αᵢ ~ 0.05 N(0, σ²ᵢ) + 0.95δ₀");
    prior_strings.emplace_back(
        fmt::format("      └─ {}", detail::sigma_prior("ᵢ", nu, s2)));
    return prior_strings;
}

void BayesBTrait::operator()(
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
    VectorXd& sigma = state.marker_variance;
    VectorXi& tracker = state.tracker;

    const auto& design_matrix = effect.design_matrix.matrix();
    const auto& cols_norm = effect.cols_norm;

    // col_norm is n - 1, since we have normalized the design matrix
    const double sqrt_sigma_e = sqrt(sigma_e);

    // Setup distributions for sampling
    std::normal_distribution<double> normal{0, 1};
    std::uniform_real_distribution<double> uniform{0, 1};
    detail::ScaledInvChiSq chi_squared{effect.prior};

    for (Index i = 0; i < coeff.size(); ++i)
    {
        if (effect.is_monomorphic(i))
        {
            continue;
        }

        // for convenience
        const double old_i = coeff(i);
        const auto& col = design_matrix.col(i);
        const double sigma_i = sigma(i);
        const double col_norm = cols_norm(i);

        const double percision_kernel = 1 / (col_norm + sigma_e / sigma_i);

        const double post_stddev = sqrt_sigma_e * sqrt(percision_kernel);
        const double logdetV = log((sigma_i * col_norm / sigma_e) + 1);

        // calculate the posterior mean and standard deviation
        double rhs = col.dot(y_adj);
        if (old_i != 0)
        {
            rhs += col_norm * old_i;
        }
        const double post_mean = rhs * percision_kernel;

        const double L_diff = (-0.5 * (logdetV - post_mean * rhs / sigma_e))
                              + logpi(1) - logpi(0);
        const double accept_prob = 1 / (1 + std::exp(L_diff));
        const int dist_index = (uniform(rng) < accept_prob) ? 0 : 1;
        tracker(i) = dist_index;

        // sample a new coefficient
        double new_i = 0.0;
        if (dist_index == 1)
        {
            new_i = (normal(rng) * post_stddev) + post_mean;
            const double diff = old_i - new_i;
            y_adj.array() += col.array() * diff;
            u.array() -= col.array() * diff;
        }
        else if (old_i != 0.0)
        {
            y_adj.array() += col.array() * old_i;
            u.array() -= col.array() * old_i;
        }
        coeff(i) = new_i;

        // Sample variance for this marker (like BayesA but with indicator)
        if (dist_index == 1)
        {
            chi_squared.compute(new_i * new_i);
            sigma(i) = chi_squared(rng);
        }
        else
        {
            // For inactive markers, keep the variance but don't update it
            // This maintains the variance structure for potential reactivation
        }
    }
    state.effect_variance = detail::var(state.u)(0);
}

}  // namespace gelex
