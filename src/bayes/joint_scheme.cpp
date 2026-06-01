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

#include <array>
#include <cmath>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/recipe_options.h"
#include "gelex/bayes/scheme.h"
#include "gelex/exception.h"
#include "gelex/types/constrained_vector.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

BayesCDScheme::BayesCDScheme(const BayesRecipeConfig& options)
    : options_{options}
{
    if (options_.modes.size() != 2)
    {
        throw GelexException("CD requires --mode AD (both effects)");
    }
    if (options_.additive.proportion() || options_.additive.proportion_update())
    {
        throw GelexException(
            "CD does not accept per-effect proportion override for A; "
            "use --joint-pi instead");
    }
    if (options_.dominance.proportion()
        || options_.dominance.proportion_update())
    {
        throw GelexException(
            "CD does not accept per-effect proportion override for D; "
            "use --joint-pi instead");
    }
    if (options_.additive.multiplier())
    {
        throw GelexException(
            "CD does not accept per-effect multiplier override for A");
    }
    if (options_.dominance.multiplier())
    {
        throw GelexException(
            "CD does not accept per-effect multiplier override for D");
    }
    if (options_.joint_proportion && options_.joint_proportion->size() != 4)
    {
        throw GelexException(
            "CD requires joint_proportion to have exactly 4 components: "
            "both-off, A-only, D-only, both-on");
    }
}

auto BayesCDScheme::make_prior(const BayesModel& model) const
    -> std::vector<GeneticPrior>
{
    const Simplex<double> proportion = options_.joint_proportion.value_or(
        Simplex<double>{{0.991, 0.003, 0.003, 0.003}});
    // proportion layout: [both-off, A-only, D-only, both-on]
    const double active_a = proportion[1] + proportion[3];
    const double active_d = proportion[2] + proportion[3];

    const double h2_a = options_.additive.heritability()
                            ? options_.additive.heritability()->value()
                            : 0.5;
    const double h2_d = options_.dominance.heritability()
                            ? options_.dominance.heritability()->value()
                            : 0.2;

    const auto* genetic_a = model.genetic(GeneticMode::A);
    if (genetic_a == nullptr)
    {
        throw GelexException("genetic design not found for mode A");
    }
    if (genetic_a->design_variance <= 0)
    {
        throw GelexException(
            "genetic design_variance must be positive for mode A");
    }
    if (active_a <= 0)
    {
        throw GelexException(
            fmt::format(
                "active_marker_weight must be positive, got {}", active_a));
    }
    const double target_a = model.phenotype_variance() * h2_a
                            / (active_a * genetic_a->design_variance);
    if (!std::isfinite(target_a) || target_a <= 0)
    {
        throw GelexException(
            fmt::format(
                "target_marker_variance must be finite and positive, got {}",
                target_a));
    }

    const auto* genetic_d = model.genetic(GeneticMode::D);
    if (genetic_d == nullptr)
    {
        throw GelexException("genetic design not found for mode D");
    }
    if (genetic_d->design_variance <= 0)
    {
        throw GelexException(
            "genetic design_variance must be positive for mode D");
    }
    if (active_d <= 0)
    {
        throw GelexException(
            fmt::format(
                "active_marker_weight must be positive, got {}", active_d));
    }
    const double target_d = model.phenotype_variance() * h2_d
                            / (active_d * genetic_d->design_variance);
    if (!std::isfinite(target_d) || target_d <= 0)
    {
        throw GelexException(
            fmt::format(
                "target_marker_variance must be finite and positive, got {}",
                target_d));
    }

    if (options_.joint_proportion_update.value_or(true))
    {
        const auto n = static_cast<Eigen::Index>(proportion.size());
        return std::vector<
            GeneticPrior>{JointGeneticPrior{JointGaussianMixturePrior{
            JointSharedMarkerVariance{std::array{
                SharedMarkerVariance{VarianceParameter{
                    target_a,
                    ScaledInvChiSqPrior{4.0, (4.0 - 2.0) / 4.0 * target_a}}},
                SharedMarkerVariance{VarianceParameter{
                    target_d,
                    ScaledInvChiSqPrior{4.0, (4.0 - 2.0) / 4.0 * target_d}}}}},
            MixtureProportion{SimplexParameter{
                proportion.to_mat(),
                DirichletPrior{Eigen::VectorXd::Ones(n)}}}}}};
    }
    return std::vector<GeneticPrior>{
        JointGeneticPrior{JointGaussianMixturePrior{
            JointSharedMarkerVariance{std::array{
                SharedMarkerVariance{VarianceParameter{
                    target_a,
                    ScaledInvChiSqPrior{4.0, (4.0 - 2.0) / 4.0 * target_a}}},
                SharedMarkerVariance{VarianceParameter{
                    target_d,
                    ScaledInvChiSqPrior{4.0, (4.0 - 2.0) / 4.0 * target_d}}}}},
            MixtureProportion{proportion.to_mat()}}}};
}

}  // namespace gelex::bayes
