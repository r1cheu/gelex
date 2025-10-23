#pragma once

#include <random>

#include "gelex/model/bayes/model.h"

namespace gelex::detail::AdditiveSampler
{

struct B
{
    auto operator()(
        const BayesModel& model,
        BayesState& states,
        std::mt19937_64& rng) const -> void;
};

}  // namespace gelex::detail::AdditiveSampler
