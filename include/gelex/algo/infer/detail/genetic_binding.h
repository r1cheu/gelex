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

#ifndef GELEX_ALGO_INFER_DETAIL_GENETIC_BINDING_H_
#define GELEX_ALGO_INFER_DETAIL_GENETIC_BINDING_H_

#include <cstddef>
#include <optional>
#include <utility>

#include <fmt/format.h>

#include "gelex/data/genotype/genotype.h"
#include "gelex/exception.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/prior_state.h"
#include "gelex/model/bayes/state.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::infer::detail
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
struct GeneticBlockDeps
{
    const bayes::GeneticEffect& effect;
    const bayes::GeneticPrior& prior;
    bayes::GeneticState& state;
    bayes::GeneticPriorState& prior_state;
    std::size_t slot;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

inline auto validate_genetic_block_shape(
    const bayes::GeneticEffect& effect,
    const bayes::GeneticState& state) -> void
{
    const auto rows = effect.X.rows();
    const auto cols = effect.X.cols();
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

struct GeneticPriorMatch
{
    const bayes::GeneticPrior* prior{};
    std::size_t block_index{};
};

inline auto find_prior_for_mode(
    const bayes::BayesPrior& prior,
    GeneticMode mode) -> std::optional<GeneticPriorMatch>
{
    std::size_t index = 0;
    for (const auto& block : prior.genetics())
    {
        if (block.contains(mode))
        {
            return GeneticPriorMatch{.prior = &block, .block_index = index};
        }
        ++index;
    }
    return std::nullopt;
}

// Binds the (effect, prior, state) trio for `mode`.
// `Ctx` must expose `model`, `prior`, and `state`.
template <typename Ctx>
inline auto bind_genetic_block(const Ctx& ctx, GeneticMode mode)
{
    const auto* effect = ctx.model.genetic(mode);
    if (effect == nullptr)
    {
        throw GelexException(
            fmt::format(
                "model has no genetic block for mode {}",
                EffectType::from_genetic(mode)));
    }
    const auto prior_match = find_prior_for_mode(ctx.prior, mode);
    if (!prior_match)
    {
        throw GelexException(
            fmt::format(
                "prior has no genetic block for mode {}",
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
    auto& block_state = ctx.state.genetic_block(prior_match->block_index);
    if (!block_state.contains(mode))
    {
        throw GelexException(
            fmt::format(
                "state has no genetic prior block for mode {}",
                EffectType::from_genetic(mode)));
    }
    validate_genetic_block_shape(*effect, *genetic_state);
    return GeneticBlockDeps{
        .effect = *effect,
        .prior = *prior_match->prior,
        .state = *genetic_state,
        .prior_state = block_state.prior_state(),
        .slot = block_state.slot(mode),
    };
}

// Binds two genetic blocks for joint updates.
template <typename Ctx>
inline auto bind_genetic_block_pair(
    const Ctx& ctx,
    GeneticMode first_mode,
    GeneticMode second_mode)
{
    auto first = bind_genetic_block(ctx, first_mode);
    auto second = bind_genetic_block(ctx, second_mode);

    if (&first.state == &second.state)
    {
        throw GelexException(
            "joint genetic sampler requires two distinct genetic blocks");
    }
    if (first.effect.X.rows() != second.effect.X.rows()
        || first.effect.X.cols() != second.effect.X.cols())
    {
        throw GelexException(
            fmt::format(
                "joint genetic blocks must have the same shape: "
                "{} is {}x{}, {} is {}x{}",
                EffectType::from_genetic(first_mode),
                first.effect.X.rows(),
                first.effect.X.cols(),
                EffectType::from_genetic(second_mode),
                second.effect.X.rows(),
                second.effect.X.cols()));
    }
    return std::pair{first, second};
}

}  // namespace gelex::infer::detail

#endif  // GELEX_ALGO_INFER_DETAIL_GENETIC_BINDING_H_
