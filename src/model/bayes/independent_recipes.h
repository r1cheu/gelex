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

#ifndef GELEX_MODEL_BAYES_INDEPENDENT_RECIPES_H_
#define GELEX_MODEL_BAYES_INDEPENDENT_RECIPES_H_

#include <memory>

#include "gelex/model/bayes/recipe_options.h"
#include "recipe_impl.h"

namespace gelex
{
class BayesModel;
}

namespace gelex::bayes
{

class BayesRRMethod final : public IndependentMethod
{
   public:
    explicit BayesRRMethod(const BayesRecipeConfig& options);

   private:
    auto make_single_genetic_prior(
        GeneticMode mode,
        const EffectConfig& effect,
        const BayesModel& model) const -> SingleGeneticPrior final;
};

class BayesAMethod final : public IndependentMethod
{
   public:
    explicit BayesAMethod(const BayesRecipeConfig& options);

   private:
    auto make_single_genetic_prior(
        GeneticMode mode,
        const EffectConfig& effect,
        const BayesModel& model) const -> SingleGeneticPrior final;
};

class BayesBMethod final : public IndependentMethod
{
   public:
    explicit BayesBMethod(const BayesRecipeConfig& options);

   private:
    auto make_single_genetic_prior(
        GeneticMode mode,
        const EffectConfig& effect,
        const BayesModel& model) const -> SingleGeneticPrior final;
};

class BayesCMethod final : public IndependentMethod
{
   public:
    explicit BayesCMethod(const BayesRecipeConfig& options);

   private:
    auto make_single_genetic_prior(
        GeneticMode mode,
        const EffectConfig& effect,
        const BayesModel& model) const -> SingleGeneticPrior final;
};

class BayesRMethod final : public IndependentMethod
{
   public:
    explicit BayesRMethod(const BayesRecipeConfig& options);

   private:
    auto make_single_genetic_prior(
        GeneticMode mode,
        const EffectConfig& effect,
        const BayesModel& model) const -> SingleGeneticPrior final;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_INDEPENDENT_RECIPES_H_
