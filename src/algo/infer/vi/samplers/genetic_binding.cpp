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

#include "gelex/algo/infer/vi/samplers/detail/genetic_binding.h"

#include <fmt/format.h>

#include "gelex/exception.h"

namespace gelex::vi::detail
{

auto bind_genetic_block(const Context& ctx, GeneticMode mode)
    -> GeneticBlockDeps
{
    const auto* effect = ctx.model.genetic(mode);
    if (effect == nullptr)
    {
        throw GelexException(
            fmt::format(
                "model has no genetic block for mode {}",
                EffectType::from_genetic(mode)));
    }
    const auto* prior = ctx.priors.genetic(mode);
    if (prior == nullptr)
    {
        throw GelexException(
            fmt::format(
                "priors has no genetic block for mode {}",
                EffectType::from_genetic(mode)));
    }
    auto* genetic_state = ctx.state.genetic(mode);
    if (genetic_state == nullptr)
    {
        throw GelexException(
            fmt::format(
                "state has no genetic block for mode {}",
                EffectType::from_genetic(mode)));
    }
    return {
        .effect = *effect,
        .prior = *prior,
        .state = *genetic_state,
    };
}

}  // namespace gelex::vi::detail
