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

#ifndef GELEX_BAYES_DETAIL_PIP_FACTORY_H_
#define GELEX_BAYES_DETAIL_PIP_FACTORY_H_

#include <cstddef>
#include <utility>

#include "gelex/bayes/genetic/detail/pip_support.h"
#include "gelex/bayes/genetic/draws.h"
#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/result.h"
#include "gelex/bayes/genetic/scaled_mixture.h"
#include "gelex/bayes/genetic/spike_slab.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/genetic_mode.h"

namespace gelex::detail
{

template <typename CoefficientDraws, typename ModeFamilyDraws>
[[nodiscard]] auto make_pip(
    const IndependentGeneticDraws<CoefficientDraws, ModeFamilyDraws>& draws)
{
    return generate_mode_values<CoefficientDraws::modes>(
        [&]<GeneticMode Mode>()
        { return make_pip(draws.template family<Mode>()); });
}

template <
    typename CoefficientDraws,
    typename ModeFamilyDraws,
    std::size_t ClassCount,
    MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_pip(
    const JointGeneticDraws<
        CoefficientDraws,
        ModeFamilyDraws,
        JointSpikeSlabDraws<ClassCount, WeightUpdate>>& draws)
{
    static_assert(ClassCount == 4);

    const auto& assignment = draws.joint_family().assignment;
    auto mode_pip = generate_mode_values<GeneticMode::A | GeneticMode::D>(
        [&]<GeneticMode Mode>()
        {
            return MarkerPipResult{assignment.probability_of(
                [](std::size_t category)
                {
                    using State = JointSpikeSlabState<ClassCount>;
                    if constexpr (Mode == GeneticMode::A)
                    {
                        return State::additive_components.at(category)
                               != State::no_component;
                    }
                    else
                    {
                        return State::dominance_components.at(category)
                               != State::no_component;
                    }
                })};
        });
    auto joint_pip
        = MarkerPipResult{assignment.probability_of(is_non_null_category)};
    return JointModeValues{std::move(mode_pip), std::move(joint_pip)};
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_DETAIL_PIP_FACTORY_H_
