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

#ifndef GELEX_ALGO_SIM_PHENOTYPE_GENERATOR_H_
#define GELEX_ALGO_SIM_PHENOTYPE_GENERATOR_H_

#include <random>

#include <Eigen/Core>

#include "gelex/algo/sim/genetic_value_calculator.h"
#include "gelex/infra/logging/simulate_event.h"

namespace gelex
{

class PhenotypeGenerator
{
   public:
    PhenotypeGenerator(
        double h2,
        double d2,
        double intercept,
        std::mt19937_64& rng);

    auto generate(
        GeneticValues& genetic_values,
        const SimulateObserver& observer = {}) -> Eigen::VectorXd;
    auto dom_scale() const -> double { return dom_scale_; }

   private:
    auto scale_dominance_if_needed(
        GeneticValues& genetic_values,
        double additive_variance,
        double d2) -> void;

    [[nodiscard]] auto calculate_residual_variance(
        double additive_variance,
        double d2) const -> double;

    auto sample_residuals(Eigen::Index n_samples, double residual_variance)
        -> Eigen::VectorXd;

    double h2_;
    double d2_;
    double intercept_;
    std::mt19937_64& rng_;
    double dom_scale_ = 1.0;
};

}  // namespace gelex

#endif  // GELEX_ALGO_SIM_PHENOTYPE_GENERATOR_H_
