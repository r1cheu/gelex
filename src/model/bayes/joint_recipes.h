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

#ifndef GELEX_MODEL_BAYES_JOINT_RECIPES_H_
#define GELEX_MODEL_BAYES_JOINT_RECIPES_H_

#include <memory>

#include "gelex/model/bayes/recipe_options.h"
#include "recipe_impl.h"

namespace gelex
{
class BayesModel;
}

namespace gelex::bayes
{

class BayesCDMethod final : public JointMethod
{
   public:
    explicit BayesCDMethod(const BayesRecipeConfig& options);

   private:
    auto make_joint_prior(
        const BayesRecipeConfig& config,
        const BayesModel& model) const -> JointGeneticPrior final;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_JOINT_RECIPES_H_
