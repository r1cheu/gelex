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

#ifndef GELEX_ALGO_INFER_VI_STEPS_RANDOM_H_
#define GELEX_ALGO_INFER_VI_STEPS_RANDOM_H_

#include <span>
#include <type_traits>

#include "gelex/algo/infer/vi/context.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/states.h"

namespace gelex::vi
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
struct RandomStepDeps
{
    std::span<const bayes::RandomEffect> effects;
    std::span<const bayes::RandomPrior> priors;
    std::span<bayes::RandomState> states;
    bayes::ResidualState& residual;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

static_assert(std::is_aggregate_v<RandomStepDeps>);

class RandomStep
{
   public:
    using Deps = RandomStepDeps;

    explicit RandomStep(Deps deps) : deps_(deps) {}

    RandomStep(const RandomStep&) = delete;
    auto operator=(const RandomStep&) -> RandomStep& = delete;
    RandomStep(RandomStep&&) noexcept = default;
    auto operator=(RandomStep&&) -> RandomStep& = delete;
    ~RandomStep() = default;

    static auto make(const Context& ctx) -> RandomStep;

    auto step() -> void;

   private:
    Deps deps_;
};

}  // namespace gelex::vi

#endif  // GELEX_ALGO_INFER_VI_STEPS_RANDOM_H_
