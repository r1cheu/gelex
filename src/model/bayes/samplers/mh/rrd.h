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

#ifndef GELEX_MODEL_BAYES_SAMPLERS_MH_RRD_H_
#define GELEX_MODEL_BAYES_SAMPLERS_MH_RRD_H_

#include <random>

namespace gelex::bayes
{
struct AdditiveEffect;
struct AdditiveState;
struct DominantEffect;
struct DominantState;
struct ResidualState;
}  // namespace gelex::bayes
namespace gelex::detail::MH
{

auto RRD(
    const bayes::AdditiveEffect& add_effect,
    bayes::AdditiveState& add_state,
    const bayes::DominantEffect& dom_effect,
    bayes::DominantState& dom_state,
    bayes::ResidualState& residual,
    std::mt19937_64& rng) -> void;
}  // namespace gelex::detail::MH

#endif  // GELEX_MODEL_BAYES_SAMPLERS_MH_RRD_H_
