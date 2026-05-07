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

#include <type_traits>
#include <utility>
#include <variant>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/algo/infer/detail/genetic_binding.h"
#include "gelex/algo/infer/mcmc/state.h"
#include "gelex/exception.h"
#include "gelex/model/bayes/model.h"

namespace gelex::mcmc
{

auto PiStep::make(const Context& ctx, GeneticMode mode) -> PiStep
{
    auto* genetic_state = ctx.state.genetic(mode);
    if (genetic_state == nullptr)
    {
        throw GelexException(
            fmt::format(
                "state has no genetic block for mode {}",
                EffectType::from_genetic(mode)));
    }
    if (!genetic_state->group.has_value())
    {
        throw GelexException(
            fmt::format(
                "PiStep requires marker allocation for genetic block {}",
                EffectType::from_genetic(mode)));
    }
    const auto* prior = infer::detail::find_prior_for_mode(ctx.method, mode);
    if (prior == nullptr)
    {
        throw GelexException(
            fmt::format(
                "method has no genetic prior for mode {}",
                EffectType::from_genetic(mode)));
    }
    if (!prior->mixture.has_value())
    {
        throw GelexException(
            fmt::format(
                "PiStep: genetic block {} has no proportion prior",
                EffectType::from_genetic(mode)));
    }
    if (!prior->mixture->proportions.estimate)
    {
        throw GelexException(
            fmt::format(
                "PiStep: genetic block {} — "
                "PiStep requires proportion.estimate = true",
                EffectType::from_genetic(mode)));
    }
    auto alpha = Eigen::VectorXd::Ones(prior->mixture->proportions.init.size());
    return PiStep{Deps{
        .group = *genetic_state->group,
        .alpha = std::move(alpha),
        .rng = ctx.rng,
    }};
}

auto PiStep::step() -> void
{
    dirichlet_.reset();

    auto update = [&](bayes::Assignment& asgn)
    { asgn.proportion = dirichlet_(asgn.count, deps_.rng); };
    std::visit(
        [&](auto& alloc)
        {
            using T = std::decay_t<decltype(alloc)>;
            if constexpr (std::is_same_v<T, bayes::Assignment>)
            {
                update(alloc);
            }
            else
            {
                update(alloc.assignment);
            }
        },
        deps_.group);
}

}  // namespace gelex::mcmc
