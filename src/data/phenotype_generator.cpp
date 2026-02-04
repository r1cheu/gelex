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

#include "gelex/data/phenotype_generator.h"

#include <algorithm>
#include <cmath>
#include <random>

#include <Eigen/Core>

#include "../src/utils/math_utils.h"

namespace gelex
{

PhenotypeGenerator::PhenotypeGenerator(PhenotypeGeneratorConfig config)
    : config_(config), rng_(config.seed)
{
}

auto PhenotypeGenerator::generate(GeneticValues& genetic_values)
    -> PhenotypeResult
{
    const double h2 = config_.h2;
    const double d2 = config_.d2;
    const double genetic_variance = detail::var(genetic_values.additive)(0);

    double residual_variance{};

    if (d2 > 0.0 && genetic_variance > 0.0)
    {
        const double dominance_raw_var
            = detail::var(genetic_values.dominance)(0);
        const double target_dominance_var = genetic_variance * d2 / h2;

        if (dominance_raw_var > 0.0)
        {
            dom_scale_ = std::sqrt(target_dominance_var / dominance_raw_var);
            genetic_values.dominance *= dom_scale_;
        }

        residual_variance = genetic_variance * (1.0 - h2 - d2) / h2;
    }
    else
    {
        residual_variance = genetic_variance * (1.0 / h2 - 1.0);
    }

    std::normal_distribution<double> residual_dist(
        0.0, std::sqrt(std::max(0.0, residual_variance)));

    Eigen::VectorXd residuals(genetic_values.additive.size());
    for (Eigen::Index i = 0; i < residuals.size(); ++i)
    {
        residuals(i) = residual_dist(rng_);
    }

    Eigen::VectorXd phenotypes
        = genetic_values.additive + genetic_values.dominance + residuals;

    if (config_.intercept != 0.0)
    {
        phenotypes.array() += config_.intercept;
    }

    double var_phen = detail::var(phenotypes)(0);

    double true_h2 = genetic_variance / var_phen;
    double true_d2 = 0.0;
    if (d2 > 0.0)
    {
        true_d2 = detail::var(genetic_values.dominance)(0) / var_phen;
    }

    return {
        .phenotypes = std::move(phenotypes),
        .true_h2 = true_h2,
        .true_d2 = true_d2};
}

}  // namespace gelex
