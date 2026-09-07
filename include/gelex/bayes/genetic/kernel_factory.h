// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_KERNEL_FACTORY_H_
#define GELEX_BAYES_GENETIC_KERNEL_FACTORY_H_

#include <cstddef>
#include <random>
#include <type_traits>
#include <utility>

#include "gelex/bayes/genetic/factory.h"
#include "gelex/bayes/genetic/gaussian_kernel.h"
#include "gelex/bayes/genetic/joint_spike_slab_kernel.h"
#include "gelex/bayes/genetic/scaled_mixture_kernel.h"
#include "gelex/bayes/genetic/spike_slab_kernel.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/state.h"
#include "gelex/genetic_mode.h"

namespace gelex
{

namespace detail
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

}  // namespace detail

template <GeneticModeSet Modes, typename... ModePriors>
[[nodiscard]] auto make_kernel(const ModeValues<Modes, ModePriors...>& prior)
{
    auto mode_kernels = transform_mode_values(
        prior,
        []<GeneticMode /*Mode*/>(const auto& mode_prior)
        { return make_kernel(mode_prior); });

    using GeneticPrior = std::remove_cvref_t<decltype(prior)>;
    using GeneticState = genetic_state_t<GeneticPrior>;
    using Sweep
        = detail::IndependentSweep<Modes, GeneticState, decltype(mode_kernels)>;
    return Sweep{std::move(mode_kernels)};
}

template <typename GeneticPrior>
using genetic_kernel_t
    = decltype(make_kernel(std::declval<const GeneticPrior&>()));

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_KERNEL_FACTORY_H_
