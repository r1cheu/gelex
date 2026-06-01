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

#ifndef GELEX_BAYES_RECIPE_H_
#define GELEX_BAYES_RECIPE_H_

#include <array>
#include <cstdint>
#include <span>
#include <string_view>
#include <utility>

#include <fmt/format.h>

#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe_options.h"
#include "gelex/bayes/scheme.h"

namespace gelex
{
class BayesModel;
}  // namespace gelex

namespace gelex::bayes
{

enum class BayesRecipeScheme : std::uint8_t
{
    RR,
    A,
    B,
    C,
    R,
    CD,
};

inline constexpr std::array BAYES_RECIPE_SCHEME_NAMES{
    std::pair{BayesRecipeScheme::RR, std::string_view{"RR"}},
    std::pair{BayesRecipeScheme::A, std::string_view{"A"}},
    std::pair{BayesRecipeScheme::B, std::string_view{"B"}},
    std::pair{BayesRecipeScheme::C, std::string_view{"C"}},
    std::pair{BayesRecipeScheme::R, std::string_view{"R"}},
    std::pair{BayesRecipeScheme::CD, std::string_view{"CD"}},
};

class BayesRecipe
{
   public:
    explicit BayesRecipe(
        BayesRecipeScheme recipe_scheme,
        BayesRecipeConfig options);
    ~BayesRecipe();

    BayesRecipe(const BayesRecipe&) = delete;
    auto operator=(const BayesRecipe&) -> BayesRecipe& = delete;

    BayesRecipe(BayesRecipe&&) noexcept = delete;
    auto operator=(BayesRecipe&&) noexcept -> BayesRecipe& = delete;

    auto make_prior(const BayesModel& model) const -> BayesPrior;

   private:
    static auto validate_modes(std::span<const GeneticMode> modes) -> void;
    auto make_random_prior(const BayesModel& model) const -> RandomPrior;
    static auto make_residual_prior(const BayesModel& model) -> ResidualPrior;

    BayesRecipeScheme recipe_scheme_;
    BayesRecipeConfig options_;
    BayesScheme scheme_;
};

auto to_bayes_recipe_scheme(std::string_view recipe_scheme)
    -> BayesRecipeScheme;

}  // namespace gelex::bayes

template <>
struct fmt::formatter<gelex::bayes::BayesRecipeScheme>
    : fmt::formatter<std::string_view>
{
    auto format(gelex::bayes::BayesRecipeScheme recipe_scheme, auto& ctx) const
    {
        return fmt::formatter<std::string_view>::format(
            to_string_view(recipe_scheme), ctx);
    }

   private:
    static constexpr auto to_string_view(
        gelex::bayes::BayesRecipeScheme recipe_scheme) -> std::string_view
    {
        for (const auto& [value, name] :
             gelex::bayes::BAYES_RECIPE_SCHEME_NAMES)
        {
            if (value == recipe_scheme)
            {
                return name;
            }
        }
        return "unknown";
    }
};

#endif  // GELEX_BAYES_RECIPE_H_
