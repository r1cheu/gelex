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
#include <string>
#include <utility>

#include <fmt/format.h>

#include "gelex/exception.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/model/bayes/prior_parameters.h"

namespace gelex::bayes
{

BayesPrior::BayesPrior(
    RandomPrior random,
    std::vector<std::unique_ptr<GeneticPrior>> genetics,
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
        auto genetic_scope = visitor.scope(GeneticPrior::name);
        for (auto [i, block] : std::views::enumerate(genetics_))
        {
            auto block_scope = visitor.scope(std::to_string(i));
            block->visit(visitor);
        }
    }
    residual_.visit(visitor);
    validate_genetics(genetics_);
}

auto BayesPrior::validate_genetics(
    const std::vector<std::unique_ptr<GeneticPrior>>& genetics) -> void
{
    std::vector<GeneticMode> seen_modes;
    for (const auto& block : genetics)
    {
        if (block == nullptr)
        {
            throw GelexException("BayesPrior: null GeneticPrior block");
        }
        const auto block_modes = block->modes();
        if (block_modes.empty())
        {
            throw GelexException("BayesPrior: GeneticPrior block has no modes");
        }
        for (const auto mode : block_modes)
        {
            if (std::ranges::contains(seen_modes, mode))
            {
                throw GelexException(
                    fmt::format(
                        "BayesPrior: duplicate GeneticMode {} across blocks",
                        mode));
            }
            seen_modes.push_back(mode);
        }
    }
}

}  // namespace gelex::bayes
