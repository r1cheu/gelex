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

#ifndef GELEX_BAYES_RECIPE_OPTIONS_H_
#define GELEX_BAYES_RECIPE_OPTIONS_H_

#include <array>
#include <cstdint>
#include <optional>
#include <string_view>
#include <utility>
#include <vector>

#include "gelex/types/constrained_value.h"
#include "gelex/types/constrained_vector.h"
#include "gelex/types/genetic_effect_type.h"

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

struct BayesRecipeOptions
{
    BayesRecipeScheme scheme{BayesRecipeScheme::RR};
    std::vector<GeneticMode> modes{GeneticMode::A};

    std::optional<OpenUnitInterval<double>> additive_heritability;
    std::optional<Simplex<double>> additive_proportion;
    std::optional<ScaleMultiplier<double>> additive_multiplier;
    std::optional<bool> additive_proportion_update;

    std::optional<OpenUnitInterval<double>> dominance_heritability;
    std::optional<Simplex<double>> dominance_proportion;
    std::optional<ScaleMultiplier<double>> dominance_multiplier;
    std::optional<bool> dominance_proportion_update;
    std::optional<OpenUnitInterval<double>> dominance_positive_probability;

    std::optional<Simplex<double>> joint_proportion;
    std::optional<bool> joint_proportion_update;

    std::optional<OpenUnitInterval<double>> random_variance_proportion;
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_RECIPE_OPTIONS_H_
