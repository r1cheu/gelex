#include "BayesRR.h"

#include <fmt/format.h>
#include <Eigen/Core>
#include <iostream>

#include "../src/logger/logger_utils.h"
#include "../src/model/bayes/bayes_effects.h"
#include "../src/model/bayes/distribution.h"
#include "data/math_utils.h"
#include "gelex/logger.h"

namespace gelex
{
using Eigen::Ref;

using Eigen::Index;

using Eigen::MatrixXd;
using Eigen::VectorXd;

VectorXd BayesRRTrait::default_sigma(Eigen::Index n_snp) const
{
    VectorXd sigma = VectorXd::Zero(1);
    return sigma;
}

VectorXd BayesRRTrait::default_pi() const
{
    VectorXd pi{{0.0, 1.0}};
    return pi;
}

bool BayesRRTrait::estimate_pi() const
{
    return false;
}

std::vector<std::string> BayesRRTrait::prior_info(
    double nu,
    double s2,
    const Ref<const VectorXd>& pi) const
{
    std::vector<std::string> prior_strings;
    prior_strings.reserve(3);
    prior_strings.emplace_back("BayesRR");
    prior_strings.emplace_back("      ├─ αᵢ ~ N(0, σ²)");
    prior_strings.emplace_back(
        fmt::format("      └─ {}", detail::sigma_prior("", nu, s2)));
    return prior_strings;
}

void BayesRRTrait::operator()(
    const bayes::AdditiveEffect& effect,
    bayes::AdditiveStatus& state,
    Ref<VectorXd> y_adj,
    double sigma_e,
    std::mt19937_64& rng) const
{
    // for convenience
    VectorXd& coeff = state.coeff;
    auto& u = state.u;
    const double old_marker_variance = state.marker_variance(0);
    const auto& design_matrix = effect.design_matrix.matrix();
    const auto& cols_norm = effect.cols_norm;

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
        const double norm = cols_norm(i);
        const double v = norm + (sigma_e / old_marker_variance);

        // calculate the posterior mean
        const double rhs = col.dot(y_adj) + (norm * old_i);
        const double post_mean = rhs / v;
        const double post_stddev = std::sqrt(sigma_e / v);

        // sample a new coefficient
        const double new_i = (normal(rng) * post_stddev) + post_mean;
        coeff(i) = new_i;

        const double diff = old_i - new_i;
        y_adj.array() += col.array() * diff;
        u.array() -= col.array() * diff;
    }

    chi_squared.compute(coeff.squaredNorm(), coeff.size() - effect.num_mono());
    state.marker_variance(0) = chi_squared(rng);
    state.effect_variance = detail::var(state.u)(0);
}

}  // namespace gelex
