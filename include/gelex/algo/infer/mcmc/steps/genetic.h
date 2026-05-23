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

#ifndef GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_H_
#define GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_H_

#include <random>
#include <type_traits>

#include "gelex/algo/infer/detail/genetic_binding.h"
#include "gelex/algo/infer/mcmc/context.h"
#include "gelex/algo/infer/mcmc/kernels/concept.h"
#include "gelex/algo/infer/mcmc/state.h"
#include "gelex/algo/infer/mcmc/sweep.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
struct GeneticStepDeps
{
    infer::detail::GeneticBlockDeps<bayes::LegacyGeneticState> block;
    bayes::ResidualState& residual;
    std::mt19937_64& rng;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

static_assert(std::is_aggregate_v<GeneticStepDeps>);

template <GeneticKernel Kernel>
class GeneticStep
{
   public:
    using Deps = GeneticStepDeps;

    explicit GeneticStep(Deps deps)
        : kernel_(deps.block.prior, deps.block.state),
          sweep_(deps.block.effect, deps.block.state, deps.residual, deps.rng)
    {
    }

    GeneticStep(const GeneticStep&) = delete;
    auto operator=(const GeneticStep&) -> GeneticStep& = delete;
    GeneticStep(GeneticStep&&) noexcept = default;
    auto operator=(GeneticStep&&) -> GeneticStep& = delete;
    ~GeneticStep() = default;

    static auto make(const Context& ctx, GeneticMode mode) -> GeneticStep
    {
        return GeneticStep{Deps{
            .block = infer::detail::bind_genetic_block(ctx, mode),
            .residual = ctx.state.residual(),
            .rng = ctx.rng,
        }};
    }

    auto step() -> void { sweep_.run(kernel_); }

   private:
    Kernel kernel_;
    GeneticSweep sweep_;
};

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_H_
