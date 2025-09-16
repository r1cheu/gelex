#include "BayesBpi.h"

#include <fmt/format.h>
#include <Eigen/Core>

#include "../src/logger/logger_utils.h"
#include "../src/model/bayes/bayes_effects.h"
#include "../src/model/bayes/distribution.h"

namespace gelex
{
using Eigen::Ref;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

bool BayesBpiTrait::estimate_pi() const
{
    return true;
}

std::vector<std::string> BayesBpiTrait::prior_info(
    double nu,
    double s2,
    const Ref<const VectorXd>& pi) const
{
    std::vector<std::string> prior_strings;
    prior_strings.reserve(4);
    prior_strings.emplace_back("BayesBπ");
    prior_strings.emplace_back("      ├─ αᵢ ~ (1-π) N(0, σ²ᵢ) + πδ₀");
    prior_strings.emplace_back(
        fmt::format("      ├─ {}", detail::sigma_prior("ᵢ", nu, s2)));
    prior_strings.emplace_back(fmt::format("      └─ π = {}", pi(0)));
    return prior_strings;
}

void BayesBpiTrait::operator()(
    const bayes::AdditiveEffect& effect,
    bayes::AdditiveStatus& state,
    Ref<VectorXd> y_adj,
    double sigma_e,
    std::mt19937_64& rng) const
{
    BayesBTrait::operator()(effect, state, y_adj, sigma_e, rng);

    state.pi.count(1) = state.tracker.sum();
    state.pi.count(0)
        = static_cast<int>(state.coeff.size() - state.pi.count(1));

    VectorXi dirichlet_counts(state.pi.count.array() + 1);
    state.pi.prop = detail::dirichlet(dirichlet_counts, rng);
}

}  // namespace gelex
