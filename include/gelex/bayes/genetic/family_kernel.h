// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_FAMILY_KERNEL_H_
#define GELEX_BAYES_GENETIC_FAMILY_KERNEL_H_

#include <type_traits>
#include <utility>

#include "gelex/bayes/genetic/construction.h"
#include "gelex/bayes/genetic/detail/independent_sweep.h"
#include "gelex/bayes/genetic/gaussian_kernel.h"
#include "gelex/bayes/genetic/joint_spike_slab_kernel.h"
#include "gelex/bayes/genetic/scaled_mixture_kernel.h"
#include "gelex/bayes/genetic/spike_slab_kernel.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/genetic_mode.h"

namespace gelex::detail
{

template <GeneticModeSet Modes, typename... ModePriors>
[[nodiscard]] auto make_kernel(const ModeValues<Modes, ModePriors...>& prior)
{
    auto mode_kernels = transform_mode_values(
        prior,
        []<GeneticMode /*Mode*/>(const auto& mode_prior)
        { return make_kernel(mode_prior); });

    using GeneticPrior = std::remove_cvref_t<decltype(prior)>;
    using GeneticState = genetic_state_t<GeneticPrior>;
    using Sweep = IndependentSweep<Modes, GeneticState, decltype(mode_kernels)>;
    return Sweep{std::move(mode_kernels)};
}

template <typename GeneticPrior>
using genetic_kernel_t
    = decltype(make_kernel(std::declval<const GeneticPrior&>()));

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_FAMILY_KERNEL_H_
