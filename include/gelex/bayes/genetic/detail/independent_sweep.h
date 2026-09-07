// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_DETAIL_INDEPENDENT_SWEEP_H_
#define GELEX_BAYES_GENETIC_DETAIL_INDEPENDENT_SWEEP_H_

#include <cstddef>
#include <random>
#include <utility>

#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/state.h"
#include "gelex/genetic_mode.h"

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

#endif  // GELEX_BAYES_GENETIC_DETAIL_INDEPENDENT_SWEEP_H_
