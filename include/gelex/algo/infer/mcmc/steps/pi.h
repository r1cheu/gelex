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

#ifndef GELEX_ALGO_INFER_MCMC_STEPS_PI_H_
#define GELEX_ALGO_INFER_MCMC_STEPS_PI_H_

#include <random>
#include <type_traits>
#include <utility>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/context.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/prior_state.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
struct PiStepDeps
{
    bayes::ProportionState& proportion;
    Eigen::VectorXd alpha;
    std::mt19937_64& rng;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

static_assert(std::is_aggregate_v<PiStepDeps>);

class PiStep
{
   public:
    using Deps = PiStepDeps;

    explicit PiStep(Deps deps) : deps_(std::move(deps)), dirichlet_(deps_.alpha)
    {
    }

    PiStep(const PiStep&) = delete;
    auto operator=(const PiStep&) -> PiStep& = delete;
    PiStep(PiStep&&) noexcept = default;
    auto operator=(PiStep&&) -> PiStep& = delete;
    ~PiStep() = default;

    static auto make(const Context& ctx, GeneticMode mode) -> PiStep;

    auto step() -> void;

   private:
    Deps deps_;
    stats::DirichletSampler<double> dirichlet_;
};

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_STEPS_PI_H_
