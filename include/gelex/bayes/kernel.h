// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_KERNEL_H_
#define GELEX_BAYES_KERNEL_H_

#include <cstddef>
#include <random>
#include <utility>
#include <vector>

#include "gelex/bayes/common_kernel.h"
#include "gelex/bayes/genetic/kernel_factory.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"

namespace gelex
{

template <typename GeneticPrior>
class BayesKernel
{
   public:
    BayesKernel(
        std::vector<RandomEffectKernel> random,
        genetic_kernel_t<GeneticPrior> genetic,
        ResidualVarianceKernel residual)
        : random_{std::move(random)},
          genetic_{std::move(genetic)},
          residual_{residual}
    {
    }

    auto step(
        const BayesModel& model,
        BayesState<GeneticPrior>& state,
        std::mt19937_64& rng) -> void
    {
        detail::update_fixed_effects(
            model.fixed(), state.fixed(), state.residual(), rng);
        for (std::size_t block = 0; block < random_.size(); ++block)
        {
            random_[block].step(
                model.random()[block],
                state.random()[block],
                state.residual(),
                rng);
        }
        genetic_.step(model.genetic(), state.genetic(), state.residual(), rng);
        residual_.step(state.residual(), rng);
    }

   private:
    std::vector<RandomEffectKernel> random_;
    genetic_kernel_t<GeneticPrior> genetic_;
    ResidualVarianceKernel residual_;
};

template <typename GeneticPrior>
[[nodiscard]] auto make_kernel(const BayesPrior<GeneticPrior>& prior)
    -> BayesKernel<GeneticPrior>
{
    std::vector<RandomEffectKernel> random;
    random.reserve(prior.random().size());
    for (const auto& parameter : prior.random())
    {
        random.emplace_back(parameter);
    }
    return BayesKernel<GeneticPrior>{
        std::move(random),
        make_kernel(prior.genetic()),
        ResidualVarianceKernel{prior.residual()}};
}

}  // namespace gelex

#endif  // GELEX_BAYES_KERNEL_H_
