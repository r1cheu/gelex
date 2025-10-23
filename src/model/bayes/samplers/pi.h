#pragma once

#include <random>

#include "gelex/model/bayes/model.h"

namespace gelex::detail::Pi
{

/**
 * @brief Sampler for pi parameter in Bayesian models
 *
 * This sampler handles the Dirichlet update for pi parameters.
 * Currently supports 2-component pi (like BayesCpi, BayesBpi),
 * but designed to be extensible for multi-component pi (like BayesR).
 */
struct Pi
{
    auto operator()(
        const BayesModel& model,
        BayesState& states,
        std::mt19937_64& rng) const -> void;
};

}  // namespace gelex::detail::Pi
