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

#ifndef GELEX_ALGO_INFER_VI_STEPS_FIXED_H_
#define GELEX_ALGO_INFER_VI_STEPS_FIXED_H_

#include <type_traits>

#include "gelex/algo/infer/vi/context.h"
#include "gelex/model/bayes/states.h"
#include "gelex/types/fixed_effects.h"

namespace gelex::vi
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
struct FixedStepDeps
{
    const FixedEffect& effect;
    bayes::FixedState& state;
    bayes::ResidualState& residual;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

static_assert(std::is_aggregate_v<FixedStepDeps>);

class FixedStep
{
   public:
    using Deps = FixedStepDeps;

    explicit FixedStep(Deps deps) : deps_(deps) {}

    FixedStep(const FixedStep&) = delete;
    auto operator=(const FixedStep&) -> FixedStep& = delete;
    FixedStep(FixedStep&&) noexcept = default;
    auto operator=(FixedStep&&) -> FixedStep& = delete;
    ~FixedStep() = default;

    static auto make(const Context& ctx) -> FixedStep
    {
        return FixedStep{Deps{
            .effect = ctx.model.fixed(),
            .state = ctx.state.fixed(),
            .residual = ctx.state.residual(),
        }};
    }

    auto step() -> void;

   private:
    Deps deps_;
};

}  // namespace gelex::vi

#endif  // GELEX_ALGO_INFER_VI_STEPS_FIXED_H_
