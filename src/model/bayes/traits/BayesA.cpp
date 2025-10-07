#include "BayesA.h"

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

VectorXd BayesATrait::default_sigma(Index n_snp) const
{
    VectorXd sigma = VectorXd::Zero(n_snp);
    return sigma;
}

VectorXd BayesATrait::default_pi() const
{
    VectorXd pi{{0.0, 1.0}};
    return pi;
}

bool BayesATrait::estimate_pi() const
{
    return false;
}

std::vector<std::string> BayesATrait::prior_info(
    double nu,
    double s2,
    const Ref<const VectorXd>& pi) const
{
    std::vector<std::string> prior_strings;
    prior_strings.reserve(3);
    prior_strings.emplace_back("BayesA");
    prior_strings.emplace_back("      ├─ αᵢ ~ N(0, σ²ᵢ)");
    prior_strings.emplace_back(
        fmt::format("      └─ {}", detail::sigma_prior("ᵢ", nu, s2)));
    return prior_strings;
}

void BayesATrait::operator()(
    const bayes::AdditiveEffect& effect,
    bayes::AdditiveStatus& state,
    Ref<VectorXd> y_adj,
    double sigma_e,
    std::mt19937_64& rng) const
{
    // for convenience
    VectorXd& coeff = state.coeff;
    auto& u = state.u;
    VectorXd& sigma = state.marker_variance;
    const auto& design_matrix = effect.design_matrix.matrix();
    const auto& cols_norm = effect.cols_norm;

    // col_norm is n - 1, since we have normalized the design matrix
    const int n = static_cast<int>(design_matrix.rows());
    const double sqrt_sigma_e = sqrt(sigma_e);

    // Setup distributions for sampling
    detail::ScaledInvChiSq chi_squared{effect.prior};
    std::normal_distribution<double> normal{0, 1};

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

        const double percision_kernel = 1 / (col_norm + sigma_e / sigma(i));

        // calculate the posterior mean and standard deviation
        const double rhs = col.dot(y_adj) + (col_norm * old_i);
        const double post_mean = rhs * percision_kernel;
        const double post_stddev = sqrt_sigma_e * sqrt(percision_kernel);

        // sample a new coefficient
        const double new_i = (normal(rng) * post_stddev) + post_mean;
        coeff(i) = new_i;

        // sample a new variance for current coefficient
        chi_squared.compute(new_i * new_i);
        sigma(i) = chi_squared(rng);

        // update the y_adj and u vectors
        const double diff = old_i - new_i;
        y_adj.array() += col.array() * diff;
        u.array() -= col.array() * diff;
    }

    state.effect_variance = detail::var(state.u)(0);
}

}  // namespace gelex
