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

#include "gelex/algo/infer/mcmc/steps/pi.h"

#include <utility>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/algo/infer/detail/genetic_binding.h"
#include "gelex/algo/infer/mcmc/state.h"
#include "gelex/exception.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/state_capabilities.h"

namespace gelex::mcmc
{

auto PiStep::make(const Context& ctx, GeneticMode mode) -> PiStep
{
    auto block = infer::detail::bind_genetic_block(ctx, mode);
    auto* cap = block.prior_state.query<bayes::ProportionStateCap>();
    if (cap == nullptr)
    {
        throw GelexException(
            fmt::format(
                "PiStep: genetic block {} has no proportion state",
                EffectType::from_genetic(mode)));
    }
    auto proportions = cap->proportion();
    if (block.slot >= proportions.size())
    {
        throw GelexException(
            fmt::format(
                "PiStep: proportion slot missing for genetic block {}",
                EffectType::from_genetic(mode)));
    }
    auto& proportion = proportions[block.slot];
    auto alpha = Eigen::VectorXd::Ones(proportion.value.size());
    return PiStep{Deps{
        .proportion = proportion,
        .alpha = std::move(alpha),
        .rng = ctx.rng,
    }};
}

auto PiStep::step() -> void
{
    dirichlet_.reset();
    if (deps_.proportion.update == bayes::ProportionUpdate::sampled)
    {
        deps_.proportion.value = dirichlet_(deps_.proportion.count, deps_.rng);
    }
}

}  // namespace gelex::mcmc
