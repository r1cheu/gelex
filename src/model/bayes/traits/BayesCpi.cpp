#include "BayesCpi.h"
#include <fmt/format.h>
#include <Eigen/Core>

#include "../src/logger/logger_utils.h"
#include "../src/model/bayes/bayes_effects.h"
#include "../src/model/bayes/distribution.h"

namespace gelex
{
using Eigen::Ref;

using Eigen::VectorXd;
using Eigen::VectorXi;

bool BayesCpiTrait::estimate_pi() const
{
    return true;
}

std::vector<std::string> BayesCpiTrait::prior_info(
    double nu,
    double s2,
    const Ref<const VectorXd>& pi) const
{
    std::vector<std::string> prior_strings;
    prior_strings.reserve(4);
    prior_strings.emplace_back("BayesCπ");
    prior_strings.emplace_back("      ├─ αᵢ ~ (1-π) N(0, σ²) + πδ₀");
    prior_strings.emplace_back(
        fmt::format("      ├─ {}", detail::sigma_prior("", nu, s2)));
    prior_strings.emplace_back(fmt::format("      └─ π = {}", pi(0)));
    return prior_strings;
}

void BayesCpiTrait::operator()(
    const bayes::AdditiveEffect& effect,
    bayes::AdditiveStatus& state,
    Ref<VectorXd> y_adj,
    double sigma_e,
    std::mt19937_64& rng) const
{
    BayesCTrait::operator()(effect, state, y_adj, sigma_e, rng);
    VectorXi dirichlet_counts(state.pi.count.array() + 1);
    state.pi.prop = detail::dirichlet(dirichlet_counts, rng);
}

}  // namespace gelex
