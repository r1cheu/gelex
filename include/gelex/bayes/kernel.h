/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef GELEX_BAYES_KERNEL_H_
#define GELEX_BAYES_KERNEL_H_

#include <cstddef>
#include <random>
#include <utility>
#include <vector>

#include "gelex/bayes/detail/common_kernel.h"
#include "gelex/bayes/detail/kernel_factory.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"

namespace gelex
{

template <typename GeneticPrior>
class BayesKernel;

template <typename GeneticPrior>
[[nodiscard]] auto make_kernel(const BayesPrior<GeneticPrior>& prior)
    -> BayesKernel<GeneticPrior>;

template <typename GeneticPrior>
class BayesKernel
{
   public:
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
    BayesKernel(
        std::vector<detail::RandomEffectKernel> random,
        detail::genetic_kernel_t<GeneticPrior> genetic,
        detail::ResidualVarianceKernel residual)
        : random_{std::move(random)},
          genetic_{std::move(genetic)},
          residual_{residual}
    {
    }

    template <typename OtherGeneticPrior>
    friend auto make_kernel(const BayesPrior<OtherGeneticPrior>& prior)
        -> BayesKernel<OtherGeneticPrior>;

    std::vector<detail::RandomEffectKernel> random_;
    detail::genetic_kernel_t<GeneticPrior> genetic_;
    detail::ResidualVarianceKernel residual_;
};

template <typename GeneticPrior>
[[nodiscard]] auto make_kernel(const BayesPrior<GeneticPrior>& prior)
    -> BayesKernel<GeneticPrior>
{
    std::vector<detail::RandomEffectKernel> random;
    random.reserve(prior.random().size());
    for (const auto& parameter : prior.random())
    {
        random.emplace_back(parameter);
    }
    return BayesKernel<GeneticPrior>{
        std::move(random),
        detail::make_kernel(prior.genetic()),
        detail::ResidualVarianceKernel{prior.residual()}};
}

}  // namespace gelex

#endif  // GELEX_BAYES_KERNEL_H_
