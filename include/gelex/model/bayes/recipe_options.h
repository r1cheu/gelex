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

#ifndef GELEX_MODEL_BAYES_RECIPE_OPTIONS_H_
#define GELEX_MODEL_BAYES_RECIPE_OPTIONS_H_

#include <optional>
#include <vector>

#include "gelex/types/constrained_value.h"
#include "gelex/types/constrained_vector.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

class EffectConfig
{
   public:
    EffectConfig() = default;
    EffectConfig(
        std::optional<OpenUnitInterval<double>> heritability,
        std::optional<Simplex<double>> proportion,
        std::optional<ScaleMultiplier<double>> multiplier,
        std::optional<bool> proportion_update);

    auto heritability() const -> const std::optional<OpenUnitInterval<double>>&
    {
        return heritability_;
    }
    auto proportion() const -> const std::optional<Simplex<double>>&
    {
        return proportion_;
    }
    auto multiplier() const -> const std::optional<ScaleMultiplier<double>>&
    {
        return multiplier_;
    }
    auto proportion_update() const -> const std::optional<bool>&
    {
        return proportion_update_;
    }

   private:
    std::optional<OpenUnitInterval<double>> heritability_;
    std::optional<Simplex<double>> proportion_;
    std::optional<ScaleMultiplier<double>> multiplier_;
    std::optional<bool> proportion_update_;
};

struct BayesRecipeConfig
{
    std::vector<GeneticMode> modes{GeneticMode::A};

    EffectConfig additive;
    EffectConfig dominance;
    std::optional<Simplex<double>> joint_proportion;
    std::optional<bool> joint_proportion_update;

    std::optional<OpenUnitInterval<double>> dominance_positive_probability;
    std::optional<OpenUnitInterval<double>> random_variance_proportion;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_RECIPE_OPTIONS_H_
