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

#ifndef GELEX_MODEL_BAYES_RECIPE_H_
#define GELEX_MODEL_BAYES_RECIPE_H_

#include <array>
#include <cstdint>
#include <memory>
#include <span>
#include <string_view>
#include <utility>

#include <fmt/format.h>

#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/recipe_options.h"

namespace gelex
{
class BayesModel;
}  // namespace gelex

namespace gelex::bayes
{

class BayesRecipeImpl;

enum class BayesRecipePreset : std::uint8_t
{
    RR,
    A,
    B,
    C,
    R,
    CD,
};

inline constexpr std::array kBayesRecipePresetNames{
    std::pair{BayesRecipePreset::RR, std::string_view{"RR"}},
    std::pair{BayesRecipePreset::A, std::string_view{"A"}},
    std::pair{BayesRecipePreset::B, std::string_view{"B"}},
    std::pair{BayesRecipePreset::C, std::string_view{"C"}},
    std::pair{BayesRecipePreset::R, std::string_view{"R"}},
    std::pair{BayesRecipePreset::CD, std::string_view{"CD"}},
};

class BayesRecipe
{
   public:
    explicit BayesRecipe(BayesRecipePreset preset, BayesRecipeConfig options);
    ~BayesRecipe();

    BayesRecipe(const BayesRecipe&) = delete;
    auto operator=(const BayesRecipe&) -> BayesRecipe& = delete;

    BayesRecipe(BayesRecipe&&) noexcept = delete;
    auto operator=(BayesRecipe&&) noexcept -> BayesRecipe& = delete;

    auto make_prior(const BayesModel& model) const -> BayesPrior;

   private:
    static auto make_impl(
        BayesRecipePreset preset,
        const BayesRecipeConfig& options) -> std::unique_ptr<BayesRecipeImpl>;
    static auto validate_modes(std::span<const GeneticMode> modes) -> void;
    auto make_random_prior(const BayesModel& model) const -> VarianceSpec;
    static auto make_residual_prior(const BayesModel& model) -> VarianceSpec;

    BayesRecipePreset preset_;
    BayesRecipeConfig options_;
    std::unique_ptr<BayesRecipeImpl> impl_;
};

auto to_bayes_recipe_preset(std::string_view preset) -> BayesRecipePreset;

}  // namespace gelex::bayes

template <>
struct fmt::formatter<gelex::bayes::BayesRecipePreset>
    : fmt::formatter<std::string_view>
{
    auto format(gelex::bayes::BayesRecipePreset preset, auto& ctx) const
    {
        return fmt::formatter<std::string_view>::format(
            to_string_view(preset), ctx);
    }

   private:
    static constexpr auto to_string_view(gelex::bayes::BayesRecipePreset preset)
        -> std::string_view
    {
        for (const auto& [value, name] : gelex::bayes::kBayesRecipePresetNames)
        {
            if (value == preset)
            {
                return name;
            }
        }
        return "unknown";
    }
};

#endif  // GELEX_MODEL_BAYES_RECIPE_H_
