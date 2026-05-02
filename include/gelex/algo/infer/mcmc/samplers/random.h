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

#ifndef GELEX_ALGO_INFER_MCMC_SAMPLERS_RANDOM_H_
#define GELEX_ALGO_INFER_MCMC_SAMPLERS_RANDOM_H_

#include <random>
#include <span>
#include <type_traits>

#include "gelex/algo/infer/mcmc/context.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/states.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
struct RandomSamplerDeps
{
    std::span<const bayes::RandomEffect> effects;
    const bayes::RandomPrior& prior;
    std::span<bayes::RandomState> states;
    bayes::ResidualState& residual;
    std::mt19937_64& rng;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

static_assert(std::is_aggregate_v<RandomSamplerDeps>);

class RandomSampler
{
   public:
    using Deps = RandomSamplerDeps;

    explicit RandomSampler(Deps deps)
        : deps_(deps),
          variance_sampler_(deps_.prior.param.nu, deps_.prior.param.s2)
    {
    }

    RandomSampler(const RandomSampler&) = delete;
    auto operator=(const RandomSampler&) -> RandomSampler& = delete;
    RandomSampler(RandomSampler&&) noexcept = default;
    auto operator=(RandomSampler&&) -> RandomSampler& = delete;
    ~RandomSampler() = default;

    static auto make(const Context& ctx) -> RandomSampler;

    auto sample() -> void;

   private:
    Deps deps_;
    NormalSampler<double> normal_{0.0};
    ScaledInvChi2Sampler<double> variance_sampler_;
};

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_SAMPLERS_RANDOM_H_
