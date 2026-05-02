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

#ifndef GELEX_ALGO_INFER_VI_SAMPLERS_RESIDUAL_H_
#define GELEX_ALGO_INFER_VI_SAMPLERS_RESIDUAL_H_

#include <type_traits>

#include "gelex/algo/infer/vi/context.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"

namespace gelex::vi
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
struct ResidualSamplerDeps
{
    const BayesModel& model;
    const bayes::ResidualPrior& prior;
    State& state;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

static_assert(std::is_aggregate_v<ResidualSamplerDeps>);

class ResidualSampler
{
   public:
    using Deps = ResidualSamplerDeps;

    explicit ResidualSampler(Deps deps)
        : deps_(deps), chi_squared_(deps_.prior.param)
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
            .model = ctx.model,
            .prior = ctx.priors.residual(),
            .state = ctx.state,
        }};
    }

    auto sample() -> void;

   private:
    Deps deps_;
    stats::detail::ScaledInvChiSq chi_squared_;
};

}  // namespace gelex::vi

#endif  // GELEX_ALGO_INFER_VI_SAMPLERS_RESIDUAL_H_
