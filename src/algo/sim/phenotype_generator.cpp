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

#include "gelex/algo/sim/phenotype_generator.h"

#include <algorithm>
#include <cmath>
#include <random>

#include <Eigen/Core>

#include "gelex/infra/utils/math_utils.h"

namespace gelex
{

PhenotypeGenerator::PhenotypeGenerator(
    double h2,
    double d2,
    double intercept,
    std::mt19937_64& rng)
    : h2_(h2), d2_(d2), intercept_(intercept), rng_(rng)
{
}

auto PhenotypeGenerator::generate(
    GeneticValues& genetic_values,
    const SimulateObserver& observer) -> Eigen::VectorXd
{
    const double additive_variance = detail::var(genetic_values.additive)(0);
    scale_dominance_if_needed(genetic_values, additive_variance, d2_);
    const double residual_variance
        = calculate_residual_variance(additive_variance, d2_);
    Eigen::VectorXd residuals
        = sample_residuals(genetic_values.additive.size(), residual_variance);

    Eigen::VectorXd phenotypes = genetic_values.additive + residuals;
    if (d2_ > 0.0)
    {
        phenotypes += genetic_values.dominance;
    }

    if (intercept_ != 0.0)
    {
        phenotypes.array() += intercept_;
    }

    const double var_phen = detail::var(phenotypes)(0);
    const double true_h2 = additive_variance / var_phen;
    const std::optional<double> true_d2
        = (d2_ > 0.0) ? std::optional<double>(
                            detail::var(genetic_values.dominance)(0) / var_phen)
                      : std::nullopt;

    if (observer)
    {
        SimulateEvent event;
        event.emplace<HeritabilityGeneratedEvent>(true_h2, true_d2);
        observer(event);
    }

    return phenotypes;
}

auto PhenotypeGenerator::scale_dominance_if_needed(
    GeneticValues& genetic_values,
    double additive_variance,
    double d2) -> void
{
    if (d2 > 0.0 && additive_variance > 0.0)
    {
        const double dominance_raw_var
            = detail::var(genetic_values.dominance)(0);
        const double target_dominance_var = additive_variance * d2 / h2_;

        if (dominance_raw_var > 0.0)
        {
            dom_scale_ = std::sqrt(target_dominance_var / dominance_raw_var);
            genetic_values.dominance *= dom_scale_;
        }
    }
}

auto PhenotypeGenerator::calculate_residual_variance(
    double additive_variance,
    double d2) const -> double
{
    return additive_variance * (1.0 - h2_ - d2) / h2_;
}

auto PhenotypeGenerator::sample_residuals(
    Eigen::Index n_samples,
    double residual_variance) -> Eigen::VectorXd
{
    if (residual_variance <= 0.0)
    {
        return Eigen::VectorXd::Zero(n_samples);
    }

    std::normal_distribution<double> residual_dist(
        0.0, std::sqrt(residual_variance));

    Eigen::VectorXd residuals(n_samples);
    for (Eigen::Index i = 0; i < residuals.size(); ++i)
    {
        residuals(i) = residual_dist(rng_);
    }
    return residuals;
}

}  // namespace gelex
