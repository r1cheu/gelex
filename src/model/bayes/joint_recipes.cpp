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

#include "joint_recipes.h"

#include <array>
#include <memory>

#include "gelex/exception.h"
#include "gelex/model/bayes/gaussian_prior.h"
#include "gelex/model/bayes/genetic_priors/gaussian.h"
#include "gelex/model/bayes/model.h"
#include "gelex/types/constrained_vector.h"

namespace gelex::bayes
{

BayesCDMethod::BayesCDMethod(const BayesRecipeConfig& options)
    : JointMethod{"CD", options}
{
    require_both_modes();
    reject_per_effect_proportion();
    reject_per_effect_multiplier();
    if (options.joint_proportion && options.joint_proportion->size() != 4)
    {
        throw GelexException(
            "CD requires joint_proportion to have exactly 4 components: "
            "both-off, A-only, D-only, both-on");
    }
}

auto BayesCDMethod::make_joint_prior(
    const BayesRecipeConfig& config,
    const BayesModel& model) const -> std::unique_ptr<GeneticPrior>
{
    const Simplex<double> proportion = config.joint_proportion.value_or(
        Simplex<double>{{0.991, 0.003, 0.003, 0.003}});
    const UpdatePolicy update = config.joint_proportion_update.value_or(true)
                                    ? UpdatePolicy::sampled
                                    : UpdatePolicy::fixed;

    // proportion layout: [both-off, A-only, D-only, both-on]
    const double active_a = proportion[1] + proportion[3];
    const double active_d = proportion[2] + proportion[3];

    const double h2_a = config.additive.heritability()
                            .value_or(default_heritability(GeneticMode::A))
                            .value();
    const double h2_d = config.dominance.heritability()
                            .value_or(default_heritability(GeneticMode::D))
                            .value();

    const double target_a = marker_variance_from_heritability(
        model, GeneticMode::A, h2_a, active_a);
    const double target_d = marker_variance_from_heritability(
        model, GeneticMode::D, h2_d, active_d);

    return std::make_unique<JointMixtureGaussianPrior>(
        std::array{GeneticMode::A, GeneticMode::D},
        std::array{
            make_marker_variance(MarkerVarianceLayout::shared, target_a),
            make_marker_variance(MarkerVarianceLayout::shared, target_d)},
        make_mixture_proportion(proportion, update));
}

auto BayesCDMethod::make_joint_prior_v2(
    const BayesRecipeConfig& config,
    const BayesModel& model) const -> std::unique_ptr<JointGeneticPrior>
{
    const Simplex<double> proportion = config.joint_proportion.value_or(
        Simplex<double>{{0.991, 0.003, 0.003, 0.003}});
    const UpdatePolicy update = config.joint_proportion_update.value_or(true)
                                    ? UpdatePolicy::sampled
                                    : UpdatePolicy::fixed;

    // proportion layout: [both-off, A-only, D-only, both-on]
    const double active_a = proportion[1] + proportion[3];
    const double active_d = proportion[2] + proportion[3];

    const double h2_a = config.additive.heritability()
                            .value_or(default_heritability(GeneticMode::A))
                            .value();
    const double h2_d = config.dominance.heritability()
                            .value_or(default_heritability(GeneticMode::D))
                            .value();

    const double target_a = marker_variance_from_heritability(
        model, GeneticMode::A, h2_a, active_a);
    const double target_d = marker_variance_from_heritability(
        model, GeneticMode::D, h2_d, active_d);

    return std::make_unique<JointGaussianMixturePrior>(
        JointSharedMarkerVariance{std::array{
            SharedMarkerVariance{VarianceParameter{
                target_a,
                ScaledInvChiSqPrior{4.0, (4.0 - 2.0) / 4.0 * target_a}}},
            SharedMarkerVariance{VarianceParameter{
                target_d,
                ScaledInvChiSqPrior{4.0, (4.0 - 2.0) / 4.0 * target_d}}}}},
        make_mixture_proportion(proportion, update));
}

}  // namespace gelex::bayes
