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

#ifndef GELEX_ALGO_INFER_MCMC_SAMPLERS_FIXED_H_
#define GELEX_ALGO_INFER_MCMC_SAMPLERS_FIXED_H_

#include <random>
#include <type_traits>

#include "gelex/algo/infer/mcmc/context.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/states.h"
#include "gelex/types/fixed_effects.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
struct FixedSamplerDeps
{
    const FixedEffect& effect;
    bayes::FixedState& state;
    bayes::ResidualState& residual;
    std::mt19937_64& rng;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

static_assert(std::is_aggregate_v<FixedSamplerDeps>);

class FixedSampler
{
   public:
    using Deps = FixedSamplerDeps;

    explicit FixedSampler(Deps deps) : deps_(deps) {}

    FixedSampler(const FixedSampler&) = delete;
    auto operator=(const FixedSampler&) -> FixedSampler& = delete;
    FixedSampler(FixedSampler&&) noexcept = default;
    auto operator=(FixedSampler&&) -> FixedSampler& = delete;
    ~FixedSampler() = default;

    static auto make(const Context& ctx) -> FixedSampler
    {
        return FixedSampler{Deps{
            .effect = ctx.model.fixed(),
            .state = ctx.state.fixed(),
            .residual = ctx.state.residual(),
            .rng = ctx.rng,
        }};
    }

    auto sample() -> void;

   private:
    Deps deps_;
    NormalSampler<double> normal_{0.0};
};

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_SAMPLERS_FIXED_H_
