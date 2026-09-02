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

#ifndef GELEX_BAYES_VARIANCE_BUDGET_H_
#define GELEX_BAYES_VARIANCE_BUDGET_H_

#include <algorithm>
#include <array>
#include <functional>
#include <utility>

#include "gelex/types/genetic_mode.h"

namespace gelex
{

// Initialization shares of phenotypic variance. Construction does not validate:
// domain checks run once in the BayesRecipe constructor together with the spec
// checks, so that a caller sees every input error at once. Until then
// residual() may be negative.
class VarianceBudget
{
   public:
    struct Shares
    {
        double additive{};
        double dominance{};
        double random{};
    };

    explicit constexpr VarianceBudget(Shares shares) noexcept
        : genetic_{shares.additive, shares.dominance}, random_{shares.random}
    {
    }

    [[nodiscard]] constexpr auto share(GeneticMode mode) const noexcept
        -> double
    {
        return genetic_[std::to_underlying(mode)];
    }

    [[nodiscard]] constexpr auto random() const noexcept -> double
    {
        return random_;
    }

    [[nodiscard]] constexpr auto residual() const noexcept -> double
    {
        return 1.0 - std::ranges::fold_left(genetic_, random_, std::plus{});
    }

   private:
    std::array<double, ALL_GENETIC_MODES.size()> genetic_;
    double random_;
};

inline constexpr double DEFAULT_ADDITIVE_SHARE = 0.5;
inline constexpr double DEFAULT_DOMINANCE_SHARE = 0.2;

// Derived from mode presence rather than a per-mode-set table, so the defaults
// cannot disagree with the mode set they are meant to serve. Yields Shares
// rather than a finished budget so that an input adapter overrides the fields
// the user gave and keeps the rest, without restating any default.
[[nodiscard]] constexpr auto default_shares(GeneticModeSet modes) noexcept
    -> VarianceBudget::Shares
{
    return {
        .additive
        = modes.contains(GeneticMode::A) ? DEFAULT_ADDITIVE_SHARE : 0.0,
        .dominance
        = modes.contains(GeneticMode::D) ? DEFAULT_DOMINANCE_SHARE : 0.0,
    };
}

}  // namespace gelex

#endif  // GELEX_BAYES_VARIANCE_BUDGET_H_
