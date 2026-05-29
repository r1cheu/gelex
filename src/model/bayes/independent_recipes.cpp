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

#include "independent_recipes.h"

#include <cstddef>
#include <memory>

#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/genetic_priors/gaussian.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/recipe_options.h"
#include "gelex/types/constrained_vector.h"
#include "gelex/types/genetic_effect_type.h"
#include "model/bayes/recipe_impl.h"

namespace gelex::bayes
{

namespace
{

auto proportion_update_from(const EffectConfig& effect) -> UpdatePolicy
{
    return effect.proportion_update().value_or(true) ? UpdatePolicy::sampled
                                                     : UpdatePolicy::fixed;
}

}  // namespace

BayesRRMethod::BayesRRMethod(const BayesRecipeConfig& options)
    : IndependentMethod{"RR", options}
{
    reject_joint_overrides();
    reject_per_effect_proportion();
    reject_per_effect_multiplier();
}

auto BayesRRMethod::make_single_genetic_prior(
    GeneticMode mode,
    const EffectConfig& effect,
    const BayesModel& model) const -> std::unique_ptr<SingleGeneticPrior>
{
    const double h2
        = effect.heritability().value_or(default_heritability(mode)).value();
    const double target
        = marker_variance_from_heritability(model, mode, h2, 1.0);
    return std::make_unique<SingleSharedGaussianPrior>(
        mode,
        SharedMarkerVariance{VarianceParameter{
            target, ScaledInvChiSqPrior{4.0, (4.0 - 2.0) / 4.0 * target}}});
}

BayesAMethod::BayesAMethod(const BayesRecipeConfig& options)
    : IndependentMethod{"A", options}
{
    reject_joint_overrides();
    reject_per_effect_proportion();
    reject_per_effect_multiplier();
}

auto BayesAMethod::make_single_genetic_prior(
    GeneticMode mode,
    const EffectConfig& effect,
    const BayesModel& model) const -> std::unique_ptr<SingleGeneticPrior>
{
    const double h2
        = effect.heritability().value_or(default_heritability(mode)).value();
    const double target
        = marker_variance_from_heritability(model, mode, h2, 1.0);
    return std::make_unique<SinglePerMarkerGaussianPrior>(
        mode,
        PerMarkerVariance{VarianceParameter{
            target, ScaledInvChiSqPrior{4.0, (4.0 - 2.0) / 4.0 * target}}});
}

BayesBMethod::BayesBMethod(const BayesRecipeConfig& options)
    : IndependentMethod{"B", options}
{
    reject_joint_overrides();
    reject_per_effect_multiplier();
}

auto BayesBMethod::make_single_genetic_prior(
    GeneticMode mode,
    const EffectConfig& effect,
    const BayesModel& model) const -> std::unique_ptr<SingleGeneticPrior>
{
    const Simplex<double> proportion
        = effect.proportion().value_or(Simplex<double>{{0.99, 0.01}});
    const double active_weight = 1.0 - proportion[0];
    const double h2
        = effect.heritability().value_or(default_heritability(mode)).value();
    const double target
        = marker_variance_from_heritability(model, mode, h2, active_weight);
    return std::make_unique<SinglePerMarkerSpikeSlabGaussianPrior>(
        mode,
        PerMarkerVariance{VarianceParameter{
            target, ScaledInvChiSqPrior{4.0, (4.0 - 2.0) / 4.0 * target}}},
        make_mixture_proportion(proportion, proportion_update_from(effect)));
}

BayesCMethod::BayesCMethod(const BayesRecipeConfig& options)
    : IndependentMethod{"C", options}
{
    reject_joint_overrides();
    reject_per_effect_multiplier();
}

auto BayesCMethod::make_single_genetic_prior(
    GeneticMode mode,
    const EffectConfig& effect,
    const BayesModel& model) const -> std::unique_ptr<SingleGeneticPrior>
{
    const Simplex<double> proportion
        = effect.proportion().value_or(Simplex<double>{{0.99, 0.01}});
    const double active_weight = 1.0 - proportion[0];
    const double h2
        = effect.heritability().value_or(default_heritability(mode)).value();
    const double target
        = marker_variance_from_heritability(model, mode, h2, active_weight);
    return std::make_unique<SingleSharedSpikeSlabGaussianPrior>(
        mode,
        SharedMarkerVariance{VarianceParameter{
            target, ScaledInvChiSqPrior{4.0, (4.0 - 2.0) / 4.0 * target}}},
        make_mixture_proportion(proportion, proportion_update_from(effect)));
}

BayesRMethod::BayesRMethod(const BayesRecipeConfig& options)
    : IndependentMethod{"R", options}
{
    reject_joint_overrides();
    require_paired_proportion_and_multiplier();
}

auto BayesRMethod::make_single_genetic_prior(
    GeneticMode mode,
    const EffectConfig& effect,
    const BayesModel& model) const -> std::unique_ptr<SingleGeneticPrior>
{
    const Simplex<double> proportion = effect.proportion().value_or(
        Simplex<double>{{0.99, 0.005, 0.003, 0.001, 0.001}});
    const ScaleMultiplier<double> multiplier = effect.multiplier().value_or(
        ScaleMultiplier<double>{{0.0, 0.001, 0.01, 0.1, 1.0}});
    double active_weight = 0.0;
    for (std::size_t i = 0; i < proportion.size(); ++i)
    {
        active_weight += proportion[i] * multiplier[i];
    }
    const double h2
        = effect.heritability().value_or(default_heritability(mode)).value();
    const double target
        = marker_variance_from_heritability(model, mode, h2, active_weight);
    return std::make_unique<SingleScaledMixtureGaussianPrior>(
        mode,
        SharedMarkerVariance{VarianceParameter{
            target, ScaledInvChiSqPrior{4.0, (4.0 - 2.0) / 4.0 * target}}},
        multiplier.to_mat(),
        make_mixture_proportion(proportion, proportion_update_from(effect)));
}

}  // namespace gelex::bayes
