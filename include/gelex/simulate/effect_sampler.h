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

#ifndef GELEX_SIMULATE_EFFECT_SAMPLER_H_
#define GELEX_SIMULATE_EFFECT_SAMPLER_H_

#include <random>
#include <span>
#include <string_view>
#include <vector>

#include "gelex/simulate/sim_types.h"

namespace gelex
{

class NormalSampler
{
   public:
    NormalSampler(
        std::span<const std::string_view> ids,
        std::span<const EffectSize> effect_sizes)
        : ids_(ids), effect_sizes_(effect_sizes)
    {
    }

    auto operator()(std::mt19937_64& rng) const -> std::vector<CausalSnp>;

   private:
    std::span<const std::string_view> ids_;
    std::span<const EffectSize> effect_sizes_;
};

}  // namespace gelex

#endif  // GELEX_SIMULATE_EFFECT_SAMPLER_H_
