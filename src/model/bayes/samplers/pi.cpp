#include "pi.h"

#include <random>

#include <Eigen/Core>

#include "gelex/model/bayes/model.h"
#include "types/bayes_effects.h"

namespace gelex::detail
{
using Eigen::VectorXd;
using Eigen::VectorXi;

namespace AdditiveSampler
{
auto Pi::operator()(
    const BayesModel& /*model*/,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    // Check if the model has additive effects with pi estimation

    if (auto* state = states.additive(); state != nullptr)
    {
        VectorXi dirichlet_counts(state->pi.count.array() + 1);

        // Sample from Dirichlet distribution
        state->pi.prop = detail::dirichlet(dirichlet_counts, rng);
    }
}

}  // namespace AdditiveSampler
//
auto DominantSampler::Pi::operator()(
    const BayesModel& /*model*/,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    // Check if the model has dominant effects with pi estimation

    if (auto* state = states.dominant(); state != nullptr)
    {
        VectorXi dirichlet_counts(state->pi.count.array() + 1);

        // Sample from Dirichlet distribution
        state->pi.prop = detail::dirichlet(dirichlet_counts, rng);
    }
}
}  // namespace gelex::detail
