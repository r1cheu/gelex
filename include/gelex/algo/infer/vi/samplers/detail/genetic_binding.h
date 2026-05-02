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

#ifndef GELEX_ALGO_INFER_VI_SAMPLERS_DETAIL_GENETIC_BINDING_H_
#define GELEX_ALGO_INFER_VI_SAMPLERS_DETAIL_GENETIC_BINDING_H_

#include <type_traits>

#include "gelex/algo/infer/vi/context.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/states.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::vi::detail
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
struct GeneticBlockDeps
{
    const bayes::GeneticEffect& effect;
    const bayes::GeneticPrior& prior;
    bayes::vi::GeneticState& state;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

static_assert(std::is_aggregate_v<GeneticBlockDeps>);

auto bind_genetic_block(const Context& ctx, GeneticMode mode)
    -> GeneticBlockDeps;

}  // namespace gelex::vi::detail

#endif  // GELEX_ALGO_INFER_VI_SAMPLERS_DETAIL_GENETIC_BINDING_H_
