#pragma once

#include <random>

#include "gelex/model/bayes/model.h"
#include "types/bayes_effects.h"
namespace gelex::detail::CommonSampler
{

struct Fixed
{
    auto operator()(
        const BayesModel& model,
        BayesState& states,
        std::mt19937_64& rng) const -> void;
};

struct Random
{
    auto operator()(
        const BayesModel& model,
        BayesState& states,
        std::mt19937_64& rng) const -> void;

   private:
    auto static sample_impl(
        const bayes::RandomEffect& effect,
        bayes::RandomState& status,
        bayes::ResidualState& residual,
        std::mt19937_64& rng) -> void;
};

struct Residual
{
    auto operator()(
        const BayesModel& model,
        BayesState& states,
        std::mt19937_64& rng) const -> void;
};

}  // namespace gelex::detail::CommonSampler
