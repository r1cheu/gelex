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

#ifndef GELEX_BAYES_SCHEME_H_
#define GELEX_BAYES_SCHEME_H_

#include <variant>
#include <vector>

#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe_options.h"

namespace gelex
{
class BayesModel;
}  // namespace gelex

namespace gelex::bayes
{

class BayesRRScheme
{
   public:
    explicit BayesRRScheme(const BayesRecipeOptions& options);
    auto make_prior(const BayesModel& model) const -> std::vector<GeneticPrior>;

   private:
    const BayesRecipeOptions& options_;
};

class BayesAScheme
{
   public:
    explicit BayesAScheme(const BayesRecipeOptions& options);
    auto make_prior(const BayesModel& model) const -> std::vector<GeneticPrior>;

   private:
    const BayesRecipeOptions& options_;
};

class BayesBScheme
{
   public:
    explicit BayesBScheme(const BayesRecipeOptions& options);
    auto make_prior(const BayesModel& model) const -> std::vector<GeneticPrior>;

   private:
    const BayesRecipeOptions& options_;
};

class BayesCScheme
{
   public:
    explicit BayesCScheme(const BayesRecipeOptions& options);
    auto make_prior(const BayesModel& model) const -> std::vector<GeneticPrior>;

   private:
    const BayesRecipeOptions& options_;
};

class BayesRScheme
{
   public:
    explicit BayesRScheme(const BayesRecipeOptions& options);
    auto make_prior(const BayesModel& model) const -> std::vector<GeneticPrior>;

   private:
    const BayesRecipeOptions& options_;
};

class BayesCDScheme
{
   public:
    explicit BayesCDScheme(const BayesRecipeOptions& options);
    auto make_prior(const BayesModel& model) const -> std::vector<GeneticPrior>;

   private:
    const BayesRecipeOptions& options_;
};

using BayesScheme = std::variant<
    BayesRRScheme,
    BayesAScheme,
    BayesBScheme,
    BayesCScheme,
    BayesRScheme,
    BayesCDScheme>;

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_SCHEME_H_
