#include "pi.h"

#include <random>

#include <Eigen/Core>

#include "../bayes_effects.h"
#include "gelex/model/bayes/model.h"

namespace gelex::detail::Pi
{
using Eigen::VectorXd;
using Eigen::VectorXi;

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

    // Get current pi counts from the state
    // For 2-component pi (like BayesCpi, BayesBpi):
    // - count(0): number of markers with coefficient = 0
    // - count(1): number of markers with coefficient != 0

    // Update counts based on current state
    // This assumes the tracker has been updated in the trait operator()
    // (as seen in BayesBpi.cpp:47-49 and BayesCpi.cpp:43-45)

    // For extensibility: this can be generalized for multi-component pi
    // by having the trait provide a count update function

    // Sample new pi proportions from Dirichlet distribution
    // Add 1 to each count for Dirichlet prior (equivalent to Beta(1,1) for 2
    // components)
}

}  // namespace gelex::detail::Pi
