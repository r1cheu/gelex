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

#ifndef GELEX_DATA_PHENOTYPE_GENERATOR_H_
#define GELEX_DATA_PHENOTYPE_GENERATOR_H_

#include <random>

#include <Eigen/Core>

#include "gelex/data/genetic_value_calculator.h"

namespace gelex
{

struct PhenotypeGeneratorConfig
{
    double h2;
    double d2;
    double intercept;
    int seed;
};

struct PhenotypeResult
{
    Eigen::VectorXd phenotypes;
    double true_h2;
    double true_d2;
};

class PhenotypeGenerator
{
   public:
    explicit PhenotypeGenerator(PhenotypeGeneratorConfig config);

    auto generate(GeneticValues& genetic_values) -> PhenotypeResult;
    auto dom_scale() const -> double { return dom_scale_; }

   private:
    PhenotypeGeneratorConfig config_;
    std::mt19937_64 rng_;
    double dom_scale_ = 1.0;
};

}  // namespace gelex

#endif  // GELEX_DATA_PHENOTYPE_GENERATOR_H_
