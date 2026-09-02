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

#ifndef GELEX_BAYES_GENETIC_KERNEL_H_
#define GELEX_BAYES_GENETIC_KERNEL_H_

#include <cstddef>
#include <random>
#include <utility>

#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/gaussian_kernel.h"
#include "gelex/bayes/genetic/independent_topology.h"
#include "gelex/bayes/genetic/joint_spike_slab_kernel.h"
#include "gelex/bayes/genetic/scaled_mixture_kernel.h"
#include "gelex/bayes/genetic/spike_slab_kernel.h"
#include "gelex/bayes/state.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::detail
{

template <GeneticModeSet Modes, typename... ModePriors>
class IndependentKernel
{
    using GeneticPrior = IndependentTopology<Modes, ModePriors...>;
    using GeneticState = genetic_state_t<GeneticPrior>;

    static auto make_mode_kernels(const GeneticPrior& prior)
    {
        return transform_mode_values(
            prior,
            []<GeneticMode /*Mode*/>(const auto& mode_prior)
            { return make_mode_kernel(mode_prior); });
    }

    using ModeKernels
        = decltype(make_mode_kernels(std::declval<const GeneticPrior&>()));

   public:
    explicit IndependentKernel(const GeneticPrior& prior)
        : mode_kernels_{make_mode_kernels(prior)}
    {
    }

    auto step(
        const bayes::GeneticDesign& design,
        GeneticState& state,
        ResidualState& residual,
        std::mt19937_64& rng) -> void
    {
        [&]<std::size_t... Index>(std::index_sequence<Index...>)
        {
            (step_mode<Modes.at(Index)>(design, state, residual, rng), ...);
        }(std::make_index_sequence<Modes.size()>{});
    }

   private:
    template <GeneticMode Mode>
    auto step_mode(
        const bayes::GeneticDesign& design,
        GeneticState& state,
        ResidualState& residual,
        std::mt19937_64& rng) -> void
    {
        mode_kernels_.template get<Mode>().template step<Mode>(
            design, state.template get<Mode>(), residual, rng);
    }

    ModeKernels mode_kernels_;
};

template <GeneticModeSet Modes, typename... ModePriors>
[[nodiscard]] auto make_genetic_kernel(
    const IndependentTopology<Modes, ModePriors...>& prior)
    -> IndependentKernel<Modes, ModePriors...>
{
    return IndependentKernel<Modes, ModePriors...>{prior};
}

template <typename GeneticPrior>
using genetic_kernel_t
    = decltype(make_genetic_kernel(std::declval<const GeneticPrior&>()));

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_KERNEL_H_
