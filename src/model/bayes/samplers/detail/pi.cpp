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

#include "gelex/model/bayes/samplers/detail/pi.h"

#include <random>
#include <type_traits>
#include <variant>

#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/samplers/detail/stick_breaking.h"
#include "gelex/model/bayes/states.h"

namespace gelex::detail
{

namespace
{

auto sample_pi(bayes::GeneticState& state, std::mt19937_64& rng) -> void
{
    if (!state.group)
    {
        return;
    }
    auto update_pi = [&](bayes::Assignment& asgn)
    {
        asgn.stick_probs = sample_stick_posteriors(asgn.count, rng);
        asgn.proportion = bayes::Assignment::stick_to_pi(asgn.stick_probs);
    };
    std::visit(
        [&](auto& alloc)
        {
            using T = std::decay_t<decltype(alloc)>;
            if constexpr (std::is_same_v<T, bayes::Assignment>)
            {
                update_pi(alloc);
            }
            else
            {
                update_pi(alloc.assignment);
            }
        },
        *state.group);
}

}  // namespace

namespace AdditiveSampler
{
auto Pi::operator()(
    const BayesModel& /*model*/,
    const bayes::Priors& /*priors*/,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    if (auto* state = states.genetic(GeneticKind::Add); state != nullptr)
    {
        sample_pi(*state, rng);
    }
}

}  // namespace AdditiveSampler

auto DominantSampler::Pi::operator()(
    const BayesModel& /*model*/,
    const bayes::Priors& /*priors*/,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    if (auto* state = states.genetic(GeneticKind::Dom); state != nullptr)
    {
        sample_pi(*state, rng);
    }
}

}  // namespace gelex::detail
