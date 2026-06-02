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

#ifndef GELEX_BAYES_MODEL_H_
#define GELEX_BAYES_MODEL_H_

#include <algorithm>
#include <vector>

#include <Eigen/Core>

#include "gelex/bayes/design.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

namespace infra
{
class FieldVisitor;
}

class BayesModel
{
   public:
    BayesModel(
        Eigen::VectorXd phenotype,
        FixedDesign fixed_design,
        std::vector<bayes::RandomDesign> random,
        std::vector<bayes::GeneticDesign> genetics);

    const FixedDesign& fixed() const { return fixed_; }

    const std::vector<bayes::RandomDesign>& random() const { return random_; }

    const std::vector<bayes::GeneticDesign>& genetics() const
    {
        return genetics_;
    }

    const bayes::GeneticDesign* genetic(GeneticMode type) const
    {
        auto it
            = std::ranges::find(genetics_, type, &bayes::GeneticDesign::type);
        return it != genetics_.end() ? &*it : nullptr;
    }

    const Eigen::VectorXd& phenotype() const { return phenotype_; }

    double phenotype_variance() const { return phenotype_var_; }
    Eigen::Index num_individuals() const { return num_individuals_; }

    auto visit(infra::FieldVisitor& visitor) const -> void;

   private:
    Eigen::Index num_individuals_{};
    double phenotype_var_{};

    Eigen::VectorXd phenotype_;

    FixedDesign fixed_;
    std::vector<bayes::RandomDesign> random_;
    std::vector<bayes::GeneticDesign> genetics_;
};

}  // namespace gelex

#endif  // GELEX_BAYES_MODEL_H_
