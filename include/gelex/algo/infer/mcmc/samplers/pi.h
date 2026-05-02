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

#ifndef GELEX_ALGO_INFER_MCMC_SAMPLERS_PI_H_
#define GELEX_ALGO_INFER_MCMC_SAMPLERS_PI_H_

#include <random>
#include <type_traits>
#include <utility>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/context.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/states.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
struct PiSamplerDeps
{
    bayes::MarkerAllocation& group;
    Eigen::VectorXd alpha;
    std::mt19937_64& rng;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

static_assert(std::is_aggregate_v<PiSamplerDeps>);

class PiSampler
{
   public:
    using Deps = PiSamplerDeps;

    explicit PiSampler(Deps deps)
        : deps_(std::move(deps)), dirichlet_(deps_.alpha)
    {
    }

    PiSampler(const PiSampler&) = delete;
    auto operator=(const PiSampler&) -> PiSampler& = delete;
    PiSampler(PiSampler&&) noexcept = default;
    auto operator=(PiSampler&&) -> PiSampler& = delete;
    ~PiSampler() = default;

    static auto make(const Context& ctx, GeneticMode mode) -> PiSampler;

    auto sample() -> void;

   private:
    Deps deps_;
    stats::DirichletSampler<double> dirichlet_;
};

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_SAMPLERS_PI_H_
