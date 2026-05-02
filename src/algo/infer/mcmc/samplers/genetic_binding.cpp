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

#include "gelex/algo/infer/mcmc/samplers/detail/genetic_binding.h"

#include <fmt/format.h>

#include "gelex/exception.h"
#include "gelex/model/bayes/genotype_storage.h"

namespace gelex::mcmc::detail
{

auto validate_genetic_block_shape(
    const bayes::GeneticEffect& effect,
    const bayes::GeneticPrior& /*prior*/,
    const bayes::GeneticState& state) -> void
{
    const auto rows = bayes::get_rows(effect.X);
    const auto cols = bayes::get_cols(effect.X);
    if (state.coeffs.size() != cols)
    {
        throw GelexException(
            fmt::format(
                "genetic block {} coeffs size {} != X cols {}",
                EffectType::from_genetic(effect.type),
                state.coeffs.size(),
                cols));
    }
    if (state.u.size() != rows)
    {
        throw GelexException(
            fmt::format(
                "genetic block {} u size {} != X rows {}",
                EffectType::from_genetic(effect.type),
                state.u.size(),
                rows));
    }
}

auto bind_genetic_block(const Context& ctx, GeneticMode mode)
    -> GeneticBlockDeps
{
    const auto* effect = ctx.model.genetic(mode);
    if (effect == nullptr)
    {
        throw GelexException(
            fmt::format(
                "model has no genetic block for mode {}",
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
    auto* genetic_state = ctx.state.genetic(mode);
    if (genetic_state == nullptr)
    {
        throw GelexException(
            fmt::format(
                "state has no genetic block for mode {}",
                EffectType::from_genetic(mode)));
    }
    validate_genetic_block_shape(*effect, *prior, *genetic_state);
    return {
        .effect = *effect,
        .prior = *prior,
        .state = *genetic_state,
    };
}

auto bind_genetic_block_pair(
    const Context& ctx,
    GeneticMode first_mode,
    GeneticMode second_mode) -> std::pair<GeneticBlockDeps, GeneticBlockDeps>
{
    auto first = bind_genetic_block(ctx, first_mode);
    auto second = bind_genetic_block(ctx, second_mode);

    if (&first.state == &second.state)
    {
        throw GelexException(
            "joint genetic sampler requires two distinct genetic blocks");
    }
    if (bayes::get_rows(first.effect.X) != bayes::get_rows(second.effect.X)
        || bayes::get_cols(first.effect.X) != bayes::get_cols(second.effect.X))
    {
        throw GelexException(
            fmt::format(
                "joint genetic blocks must have the same shape: "
                "{} is {}x{}, {} is {}x{}",
                EffectType::from_genetic(first_mode),
                bayes::get_rows(first.effect.X),
                bayes::get_cols(first.effect.X),
                EffectType::from_genetic(second_mode),
                bayes::get_rows(second.effect.X),
                bayes::get_cols(second.effect.X)));
    }
    return {first, second};
}

}  // namespace gelex::mcmc::detail
