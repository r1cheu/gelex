#pragma once

#include <random>

#include "gelex/model/bayes/model.h"

namespace gelex::detail::AdditiveSampler
{

struct RRD
{
    auto operator()(
        const BayesModel& model,
        BayesState& states,
        std::mt19937_64& rng) const -> void;

    static auto g_ad(double w_j, double a, double d) -> double
    {
        if (std::abs(w_j) < 1e-12)
        {
            return 0.5;
        }
        return (1 - w_j * sign(a) * sign(d)) / 2;
    };

    static double sign(double x) { return (x > 0.0) ? 1.0 : -1.0; }
};

}  // namespace gelex::detail::AdditiveSampler
