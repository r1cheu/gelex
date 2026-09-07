// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_VARIANCE_BUDGET_H_
#define GELEX_BAYES_VARIANCE_BUDGET_H_

#include <algorithm>
#include <array>
#include <functional>
#include <utility>

#include "gelex/genetic_mode.h"

namespace gelex
{

/**
 * @brief Phenotypic variance proportions allocated to model components.
 *
 * All supplied proportions are finite and non-negative, with
 * @f[ p_e = 1 - \sum_{m \in \mathcal M} p_m - p_r > 0. @f]
 */
class VarianceBudget
{
   public:
    struct Proportion
    {
        double additive{};
        double dominance{};
        double random{};
    };

    /**
     * @throws GelexException if a proportion is non-finite or negative, or
     * their sum is not less than one.
     */
    explicit VarianceBudget(Proportion shares);

    constexpr auto genetic(GeneticMode mode) const noexcept -> double
    {
        return genetic_[std::to_underlying(mode)];
    }

    constexpr auto random() const noexcept -> double { return random_; }

    constexpr auto residual() const noexcept -> double
    {
        return 1.0 - std::ranges::fold_left(genetic_, random_, std::plus{});
    }

   private:
    std::array<double, all_genetic_modes.size()> genetic_;
    double random_;
};

inline constexpr double default_additive_share = 0.5;
inline constexpr double default_dominance_share = 0.2;

[[nodiscard]] constexpr auto default_proportion(GeneticModeSet modes) noexcept
    -> VarianceBudget::Proportion
{
    return {
        .additive
        = modes.contains(GeneticMode::A) ? default_additive_share : 0.0,
        .dominance
        = modes.contains(GeneticMode::D) ? default_dominance_share : 0.0,
    };
}

}  // namespace gelex

#endif  // GELEX_BAYES_VARIANCE_BUDGET_H_
