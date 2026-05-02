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

#ifndef GELEX_ALGO_INFER_MCMC_SAMPLERS_RESIDUAL_H_
#define GELEX_ALGO_INFER_MCMC_SAMPLERS_RESIDUAL_H_

#include <random>
#include <type_traits>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/context.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/states.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
struct ResidualSamplerDeps
{
    Eigen::Index num_individuals;
    const bayes::ResidualPrior& prior;
    bayes::ResidualState& state;
    std::mt19937_64& rng;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

static_assert(std::is_aggregate_v<ResidualSamplerDeps>);

class ResidualSampler
{
   public:
    using Deps = ResidualSamplerDeps;

    explicit ResidualSampler(Deps deps)
        : deps_(deps),
          variance_sampler_(deps_.prior.param.nu, deps_.prior.param.s2)
    {
    }

    ResidualSampler(const ResidualSampler&) = delete;
    auto operator=(const ResidualSampler&) -> ResidualSampler& = delete;
    ResidualSampler(ResidualSampler&&) noexcept = default;
    auto operator=(ResidualSampler&&) -> ResidualSampler& = delete;
    ~ResidualSampler() = default;

    static auto make(const Context& ctx) -> ResidualSampler
    {
        return ResidualSampler{Deps{
            .num_individuals = ctx.model.num_individuals(),
            .prior = ctx.priors.residual(),
            .state = ctx.state.residual(),
            .rng = ctx.rng,
        }};
    }

    auto sample() -> void;

   private:
    Deps deps_;
    stats::ScaledInvChi2Sampler<double> variance_sampler_;
};

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_SAMPLERS_RESIDUAL_H_
