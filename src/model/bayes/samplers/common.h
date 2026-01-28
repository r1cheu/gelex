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

#ifndef GELEX_MODEL_BAYES_SAMPLERS_COMMON_H_
#define GELEX_MODEL_BAYES_SAMPLERS_COMMON_H_

#include <random>

#include "../src/types/bayes_effects.h"
#include "gelex/model/bayes/model.h"
namespace gelex::detail::CommonSampler
{

struct Fixed
{
    auto operator()(
        const BayesModel& model,
        BayesState& states,
        std::mt19937_64& rng) const -> void;
};

struct Random
{
    auto operator()(
        const BayesModel& model,
        BayesState& states,
        std::mt19937_64& rng) const -> void;

   private:
    auto static sample_impl(
        const bayes::RandomEffect& effect,
        bayes::RandomState& status,
        bayes::ResidualState& residual,
        std::mt19937_64& rng) -> void;
};

struct Residual
{
    auto operator()(
        const BayesModel& model,
        BayesState& states,
        std::mt19937_64& rng) const -> void;
};

}  // namespace gelex::detail::CommonSampler

#endif  // GELEX_MODEL_BAYES_SAMPLERS_COMMON_H_
