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

#ifndef GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_JOINT_H_
#define GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_JOINT_H_

#include <random>
#include <type_traits>

#include "gelex/algo/infer/detail/genetic_binding.h"
#include "gelex/algo/infer/mcmc/context.h"
#include "gelex/algo/infer/mcmc/joint_sweep.h"
#include "gelex/algo/infer/mcmc/kernels/concept.h"
#include "gelex/algo/infer/mcmc/state.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
struct GeneticJointStepDeps
{
    infer::detail::GeneticBlockDeps<bayes::LegacyGeneticState> first;
    infer::detail::GeneticBlockDeps<bayes::LegacyGeneticState> second;
    bayes::ResidualState& residual;
    std::mt19937_64& rng;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

static_assert(std::is_aggregate_v<GeneticJointStepDeps>);

template <GeneticJointKernel Kernel>
class GeneticJointStep
{
   public:
    using Deps = GeneticJointStepDeps;

    explicit GeneticJointStep(Deps deps)
        : kernel_(deps.first, deps.second),
          sweep_(deps.first, deps.second, deps.residual, deps.rng)
    {
    }

    GeneticJointStep(const GeneticJointStep&) = delete;
    auto operator=(const GeneticJointStep&) -> GeneticJointStep& = delete;
    GeneticJointStep(GeneticJointStep&&) noexcept = default;
    auto operator=(GeneticJointStep&&) -> GeneticJointStep& = delete;
    ~GeneticJointStep() = default;

    static auto make(
        const Context& ctx,
        GeneticMode first_mode,
        GeneticMode second_mode) -> GeneticJointStep
    {
        auto [first, second] = infer::detail::bind_genetic_block_pair(
            ctx, first_mode, second_mode);
        return GeneticJointStep{Deps{
            .first = first,
            .second = second,
            .residual = ctx.state.residual(),
            .rng = ctx.rng,
        }};
    }

    auto step() -> void { sweep_.run(kernel_); }

   private:
    Kernel kernel_;
    GeneticJointSweep sweep_;
};

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_JOINT_H_
