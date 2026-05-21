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

#ifndef GELEX_ALGO_INFER_VI_STEPS_RESIDUAL_H_
#define GELEX_ALGO_INFER_VI_STEPS_RESIDUAL_H_

#include <type_traits>

#include "gelex/algo/infer/vi/context.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/legacy_prior.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"

namespace gelex::vi
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
struct ResidualStepDeps
{
    const BayesModel& model;
    const bayes::OldVarianceSpec& variance;
    State& state;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

static_assert(std::is_aggregate_v<ResidualStepDeps>);

class ResidualStep
{
   public:
    using Deps = ResidualStepDeps;

    explicit ResidualStep(Deps deps)
        : deps_(deps),
          chi_squared_(
              stats::detail::ScaledInvChiSqParams{
                  .nu = deps_.variance.prior.degrees_of_freedom(),
                  .s2 = deps_.variance.prior.scale()})
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
            .model = ctx.model,
            .variance = ctx.method.residual,
            .state = ctx.state,
        }};
    }

    auto step() -> void;

   private:
    Deps deps_;
    stats::detail::ScaledInvChiSq chi_squared_;
};

}  // namespace gelex::vi

#endif  // GELEX_ALGO_INFER_VI_STEPS_RESIDUAL_H_
