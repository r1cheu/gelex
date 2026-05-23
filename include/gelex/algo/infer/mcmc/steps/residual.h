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

#ifndef GELEX_ALGO_INFER_MCMC_STEPS_RESIDUAL_H_
#define GELEX_ALGO_INFER_MCMC_STEPS_RESIDUAL_H_

#include <random>
#include <type_traits>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/context.h"
#include "gelex/algo/infer/mcmc/state.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/prior.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
struct ResidualStepDeps
{
    Eigen::Index num_individuals;
    const bayes::VarianceSpec& variance;
    bayes::ResidualState& state;
    std::mt19937_64& rng;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

static_assert(std::is_aggregate_v<ResidualStepDeps>);

class ResidualStep
{
   public:
    using Deps = ResidualStepDeps;

    explicit ResidualStep(Deps deps)
        : deps_(deps),
          variance_sampler_(
              deps_.variance.prior().degrees_of_freedom(),
              deps_.variance.prior().scale())
    {
    }

    ResidualStep(const ResidualStep&) = delete;
    auto operator=(const ResidualStep&) -> ResidualStep& = delete;
    ResidualStep(ResidualStep&&) noexcept = default;
    auto operator=(ResidualStep&&) -> ResidualStep& = delete;
    ~ResidualStep() = default;

    static auto make(const Context& ctx) -> ResidualStep
    {
        return ResidualStep{Deps{
            .num_individuals = ctx.model.num_individuals(),
            .variance = ctx.prior.residual(),
            .state = ctx.state.residual(),
            .rng = ctx.rng,
        }};
    }

    auto step() -> void;

   private:
    Deps deps_;
    stats::ScaledInvChi2Sampler<double> variance_sampler_;
};

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_STEPS_RESIDUAL_H_
