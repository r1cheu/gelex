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

#ifndef GELEX_MODEL_BAYES_MODEL_H_
#define GELEX_MODEL_BAYES_MODEL_H_

#include <algorithm>
#include <vector>

#include <Eigen/Core>

#include "gelex/model/bayes/effects.h"
#include "gelex/types/fixed_effects.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

class BayesModel
{
   public:
    BayesModel(
        Eigen::VectorXd phenotype,
        FixedEffect fixed_effects,
        std::vector<bayes::RandomEffect> random,
        std::vector<bayes::GeneticEffect> genetics);

    const FixedEffect& fixed() const { return fixed_; }

    const std::vector<bayes::RandomEffect>& random() const { return random_; }

    const std::vector<bayes::GeneticEffect>& genetics() const
    {
        return genetics_;
    }

    const bayes::GeneticEffect* genetic(GeneticMode type) const
    {
        auto it
            = std::ranges::find(genetics_, type, &bayes::GeneticEffect::type);
        return it != genetics_.end() ? &*it : nullptr;
    }

    const Eigen::VectorXd& phenotype() const { return phenotype_; }

    double phenotype_variance() const { return phenotype_var_; }
    Eigen::Index num_individuals() const { return num_individuals_; }

   private:
    Eigen::Index num_individuals_{};
    double phenotype_var_{};

    Eigen::VectorXd phenotype_;

    FixedEffect fixed_;
    std::vector<bayes::RandomEffect> random_;
    std::vector<bayes::GeneticEffect> genetics_;
};

}  // namespace gelex

#endif  // GELEX_MODEL_BAYES_MODEL_H_
