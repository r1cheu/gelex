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
#include <vector>

#include "gelex/bayes/common_kernel.h"
#include "gelex/bayes/genetic/kernel.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"
#include "gelex/exception.h"
#include "gelex/types/genetic_mode.h"

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
        validate_random_effects(model, state);
        validate_design(model.genetic());
        fixed_.step(model.fixed(), state.fixed(), state.residual(), rng);
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
    auto validate_random_effects(
        const BayesModel& model,
        const BayesState<GeneticPrior>& state) const -> void
    {
        if (model.random().size() != random_.size()
            || state.random().size() != random_.size())
        {
            throw GelexException(
                "random prior, design and state counts must match");
        }
    }

    static auto validate_design(const bayes::GeneticDesign& design) -> void
    {
        for (const GeneticMode mode : GeneticPrior::modes.each())
        {
            if (!design.contains(mode))
            {
                throw GelexException(
                    "genetic design does not contain every mode required by "
                    "the kernel");
            }
        }
    }

    explicit BayesKernel(const BayesPrior<GeneticPrior>& prior)
        : genetic_{detail::make_genetic_kernel(prior.genetic())},
          residual_{prior.residual()}
    {
        random_.reserve(prior.random().size());
        for (const auto& parameter : prior.random())
        {
            random_.emplace_back(parameter);
        }
    }

    template <typename T>
    friend auto make_kernel(const BayesPrior<T>& prior) -> BayesKernel<T>;

    detail::FixedEffectKernel fixed_;
    std::vector<detail::RandomEffectKernel> random_;
    detail::genetic_kernel_t<GeneticPrior> genetic_;
    detail::ResidualVarianceKernel residual_;
};

template <typename GeneticPrior>
[[nodiscard]] auto make_kernel(const BayesPrior<GeneticPrior>& prior)
    -> BayesKernel<GeneticPrior>
{
    return BayesKernel<GeneticPrior>{prior};
}

}  // namespace gelex

#endif  // GELEX_BAYES_KERNEL_H_
