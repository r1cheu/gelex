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

#include "gelex/model/bayes/prior.h"

#include <algorithm>
#include <ranges>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

#include <fmt/format.h>

#include "gelex/exception.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/model/bayes/prior_parameters.h"

namespace gelex::bayes
{

BayesPrior::BayesPrior(
    RandomPrior random,
    std::vector<GeneticPriorBlock> genetics,
    ResidualPrior residual)
    : random_(random), genetics_(std::move(genetics)), residual_(residual)
{
    validate_genetics(genetics_);
}

auto BayesPrior::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    random_.visit(visitor);
    {
        auto genetic_scope = visitor.scope("genetic");
        for (auto [i, block] : std::views::enumerate(genetics_))
        {
            auto block_scope = visitor.scope(std::to_string(i));
            std::visit(
                [&visitor](auto& prior) { prior->visit(visitor); }, block);
        }
    }
    residual_.visit(visitor);
    // FieldVisitor may mutate prior fields, so restore invariant checks here.
    validate_genetics(genetics_);
}

auto BayesPrior::validate_genetics(
    const std::vector<GeneticPriorBlock>& genetics) -> void
{
    std::set<GeneticMode> seen_modes;
    auto add_mode = [&seen_modes](GeneticMode mode)
    {
        if (!seen_modes.insert(mode).second)
        {
            throw GelexException(
                fmt::format(
                    "BayesPrior: duplicate GeneticMode {} across blocks",
                    mode));
        }
    };
    for (const auto& block : genetics)
    {
        std::visit(
            [&add_mode](const auto& prior)
            {
                using Prior = std::decay_t<decltype(prior)>;
                if (prior == nullptr)
                {
                    throw GelexException("BayesPrior: null genetic block");
                }

                if constexpr (
                    std::is_same_v<Prior, std::unique_ptr<SingleGeneticPrior>>)
                {
                    add_mode(prior->mode());
                }
                else if constexpr (
                    std::is_same_v<Prior, std::unique_ptr<JointGeneticPrior>>)
                {
                    // A joint prior owns both genetic modes, so reserve both
                    // slots for duplicate checks.
                    add_mode(GeneticMode::A);
                    add_mode(GeneticMode::D);
                }
            },
            block);
    }
}

}  // namespace gelex::bayes
