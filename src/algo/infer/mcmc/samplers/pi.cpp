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

#include "gelex/algo/infer/mcmc/samplers/pi.h"

#include <type_traits>
#include <utility>
#include <variant>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/states.h"

namespace gelex::mcmc
{

auto PiSampler::make(const Context& ctx, GeneticMode mode) -> PiSampler
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
                "PiSampler requires marker allocation for genetic block {}",
                EffectType::from_genetic(mode)));
    }
    const auto* prior = ctx.priors.genetic(mode);
    if (prior == nullptr)
    {
        throw GelexException(
            fmt::format(
                "priors has no genetic block for mode {}",
                EffectType::from_genetic(mode)));
    }
    auto alpha = std::visit(
        [&](const auto& marker_prior) -> Eigen::VectorXd
        {
            using T = std::decay_t<decltype(marker_prior)>;
            if constexpr (
                std::is_same_v<T, bayes::SpikePrior>
                || std::is_same_v<T, bayes::MixturePrior>)
            {
                if (!marker_prior.proportion.estimate)
                {
                    throw GelexException(
                        fmt::format(
                            "PiSampler: genetic block {} — "
                            "PiSampler requires proportion.estimate = true",
                            EffectType::from_genetic(mode)));
                }
                return Eigen::VectorXd::Ones(
                    marker_prior.proportion.init.size());
            }
            else
            {
                throw GelexException(
                    fmt::format(
                        "PiSampler: genetic block {} has no proportion prior",
                        EffectType::from_genetic(mode)));
            }
        },
        prior->marker);
    return PiSampler{Deps{
        .group = *genetic_state->group,
        .alpha = std::move(alpha),
        .rng = ctx.rng,
    }};
}

auto PiSampler::sample() -> void
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
