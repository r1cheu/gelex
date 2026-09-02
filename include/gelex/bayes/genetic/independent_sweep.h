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

#ifndef GELEX_BAYES_GENETIC_INDEPENDENT_SWEEP_H_
#define GELEX_BAYES_GENETIC_INDEPENDENT_SWEEP_H_

#include <cstddef>
#include <random>
#include <utility>

#include "gelex/bayes/design_data.h"
#include "gelex/bayes/state.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::detail
{

template <GeneticModeSet Modes, typename GeneticState, typename ModeKernels>
class IndependentSweep
{
   public:
    explicit IndependentSweep(ModeKernels mode_kernels)
        : mode_kernels_{std::move(mode_kernels)}
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

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_INDEPENDENT_SWEEP_H_
