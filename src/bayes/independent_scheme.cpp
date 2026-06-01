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

#include <vector>

#include "bayes/detail/scheme_helpers.h"
#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/recipe_options.h"
#include "gelex/bayes/scheme.h"
#include "gelex/types/constrained_vector.h"

namespace gelex::bayes
{
BayesRRScheme::BayesRRScheme(const BayesRecipeConfig& options)
    : options_{options}
{
    detail::reject_joint_options(options_, "RR");
    detail::reject_proportion_options(options_, "RR");
    detail::reject_multiplier_options(options_, "RR");
}

auto BayesRRScheme::make_prior(const BayesModel& model) const
    -> std::vector<GeneticPrior>
{
    std::vector<GeneticPrior> priors;
    for (const auto mode : options_.modes)
    {
        const auto& effect = detail::effect_config(options_, mode);
        const double target = detail::target_marker_variance(
            model, mode, detail::heritability(mode, effect), 1.0);

        priors.emplace_back(
            SingleGeneticPrior{SingleSharedGaussianPrior{
                mode,
                SharedMarkerVariance{detail::variance_parameter(target)}}});
    }
    return priors;
}

BayesAScheme::BayesAScheme(const BayesRecipeConfig& options) : options_{options}
{
    detail::reject_joint_options(options_, "A");
    detail::reject_proportion_options(options_, "A");
    detail::reject_multiplier_options(options_, "A");
}

auto BayesAScheme::make_prior(const BayesModel& model) const
    -> std::vector<GeneticPrior>
{
    std::vector<GeneticPrior> priors;
    for (const auto mode : options_.modes)
    {
        const auto& effect = detail::effect_config(options_, mode);
        const double target = detail::target_marker_variance(
            model, mode, detail::heritability(mode, effect), 1.0);

        priors.emplace_back(
            SingleGeneticPrior{SinglePerMarkerGaussianPrior{
                mode, PerMarkerVariance{detail::variance_parameter(target)}}});
    }
    return priors;
}

BayesBScheme::BayesBScheme(const BayesRecipeConfig& options) : options_{options}
{
    detail::reject_joint_options(options_, "B");
    detail::reject_multiplier_options(options_, "B");
}

auto BayesBScheme::make_prior(const BayesModel& model) const
    -> std::vector<GeneticPrior>
{
    std::vector<GeneticPrior> priors;
    for (const auto mode : options_.modes)
    {
        const auto& effect = detail::effect_config(options_, mode);
        const Simplex<double> proportion
            = effect.proportion().value_or(Simplex<double>{{0.99, 0.01}});
        const double target = detail::target_marker_variance(
            model,
            mode,
            detail::heritability(mode, effect),
            1.0 - proportion[0]);

        priors.emplace_back(
            SingleGeneticPrior{SinglePerMarkerSpikeSlabGaussianPrior{
                mode,
                PerMarkerVariance{detail::variance_parameter(target)},
                detail::mixture_proportion(
                    proportion, effect.proportion_update().value_or(true))}});
    }
    return priors;
}

BayesCScheme::BayesCScheme(const BayesRecipeConfig& options) : options_{options}
{
    detail::reject_joint_options(options_, "C");
    detail::reject_multiplier_options(options_, "C");
}

auto BayesCScheme::make_prior(const BayesModel& model) const
    -> std::vector<GeneticPrior>
{
    std::vector<GeneticPrior> priors;
    for (const auto mode : options_.modes)
    {
        const auto& effect = detail::effect_config(options_, mode);
        const Simplex<double> proportion
            = effect.proportion().value_or(Simplex<double>{{0.99, 0.01}});
        const double target = detail::target_marker_variance(
            model,
            mode,
            detail::heritability(mode, effect),
            1.0 - proportion[0]);

        priors.emplace_back(
            SingleGeneticPrior{SingleSharedSpikeSlabGaussianPrior{
                mode,
                SharedMarkerVariance{detail::variance_parameter(target)},
                detail::mixture_proportion(
                    proportion, effect.proportion_update().value_or(true))}});
    }
    return priors;
}

BayesRScheme::BayesRScheme(const BayesRecipeConfig& options) : options_{options}
{
    detail::reject_joint_options(options_, "R");
    detail::reject_unpaired_proportion_multiplier(options_, "R");
}

auto BayesRScheme::make_prior(const BayesModel& model) const
    -> std::vector<GeneticPrior>
{
    std::vector<GeneticPrior> priors;
    for (const auto mode : options_.modes)
    {
        const auto& effect = detail::effect_config(options_, mode);
        const Simplex<double> proportion = effect.proportion().value_or(
            Simplex<double>{{0.99, 0.005, 0.003, 0.001, 0.001}});
        const ScaleMultiplier<double> multiplier = effect.multiplier().value_or(
            ScaleMultiplier<double>{{0.0, 0.001, 0.01, 0.1, 1.0}});
        const double target = detail::target_marker_variance(
            model,
            mode,
            detail::heritability(mode, effect),
            detail::scaled_active_marker_weight(proportion, multiplier));

        priors.emplace_back(
            SingleGeneticPrior{SingleScaledMixtureGaussianPrior{
                mode,
                SharedMarkerVariance{detail::variance_parameter(target)},
                multiplier.to_mat(),
                detail::mixture_proportion(
                    proportion, effect.proportion_update().value_or(true))}});
    }
    return priors;
}

}  // namespace gelex::bayes
