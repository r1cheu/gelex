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

#ifndef GELEX_BAYES_DETAIL_KERNEL_FACTORY_H_
#define GELEX_BAYES_DETAIL_KERNEL_FACTORY_H_

#include <cstddef>
#include <type_traits>
#include <utility>

#include "gelex/bayes/detail/state_factory.h"
#include "gelex/bayes/genetic/kernel/gaussian.h"
#include "gelex/bayes/genetic/kernel/independent_sweep.h"
#include "gelex/bayes/genetic/kernel/joint_spike_slab.h"
#include "gelex/bayes/genetic/kernel/scaled_mixture.h"
#include "gelex/bayes/genetic/kernel/spike_slab.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/genetic_mode.h"

namespace gelex::detail
{

template <VarianceLayout Kind>
[[nodiscard]] auto make_kernel(const GaussianPrior<Kind>& prior)
{
    return GaussianKernel<Kind>{prior};
}

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_kernel(const SpikeSlabPrior<Kind, WeightUpdate>& prior)
{
    return SpikeSlabKernel<Kind, WeightUpdate>{prior};
}

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_kernel(
    const ScaledMixturePrior<ClassCount, WeightUpdate>& prior)
{
    return ScaledMixtureKernel<ClassCount, WeightUpdate>{prior};
}

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

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_kernel(
    const JointModeValues<
        ModeValues<
            GeneticMode::A | GeneticMode::D,
            GaussianPrior<VarianceLayout::Pooled>,
            HalfNormalPrior>,
        JointSpikeSlabPrior<ClassCount, WeightUpdate>>& prior)
{
    return JointSpikeSlabKernel<ClassCount, WeightUpdate>{prior};
}

template <typename GeneticPrior>
using genetic_kernel_t
    = decltype(make_kernel(std::declval<const GeneticPrior&>()));

}  // namespace gelex::detail

#endif  // GELEX_BAYES_DETAIL_KERNEL_FACTORY_H_
