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

#ifndef GELEX_MODEL_BAYES_PRIOR_CONFIG_H_
#define GELEX_MODEL_BAYES_PRIOR_CONFIG_H_

#include <span>
#include <vector>

#include <Eigen/Core>

#include "gelex/types/bayes_method.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

struct GeneticPriorConfig
{
    GeneticKind type;
    double heritability;
    Eigen::VectorXd proportion;
    Eigen::VectorXd multiplier;
};

namespace bayes
{
struct PriorSet;
}

class PriorSetConfig
{
   public:
    PriorSetConfig(BayesMethodConfig method, double phenotype_variance);

    auto override_proportion(GeneticKind type, std::span<const double> values)
        -> PriorSetConfig&;
    auto override_multiplier(GeneticKind type, std::span<const double> values)
        -> PriorSetConfig&;
    auto override_positive_prob(double value) -> PriorSetConfig&;

   private:
    friend struct bayes::PriorSet;

    [[nodiscard]] auto find_genetic(GeneticKind type) -> GeneticPriorConfig*;
    [[nodiscard]] auto find_genetic(GeneticKind type) const
        -> const GeneticPriorConfig*;

    BayesMethodConfig method_;
    double phenotype_variance_;
    std::vector<GeneticPriorConfig> genetics_;
    double random_variance_proportion_{0.1};
    double residual_variance_proportion_{0.3};
    double positive_prob_{0.5};
};

}  // namespace gelex

#endif  // GELEX_MODEL_BAYES_PRIOR_CONFIG_H_
