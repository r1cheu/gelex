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

#include <cmath>
#include <cstddef>
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

BayesRRScheme::BayesRRScheme(const BayesRecipeConfig& options)
    : options_{options}
{
    if (options_.joint_proportion)
    {
        throw GelexException("RR does not accept --joint-pi");
    }
    if (options_.joint_proportion_update)
    {
        throw GelexException("RR does not accept --estimate-joint-pi");
    }
    for (const auto mode : options_.modes)
    {
        const auto& effect
            = mode == GeneticMode::A ? options_.additive : options_.dominance;
        if (effect.proportion() || effect.proportion_update())
        {
            throw GelexException(
                fmt::format(
                    "RR does not accept proportion override for {}", mode));
        }
        if (effect.multiplier())
        {
            throw GelexException(
                fmt::format(
                    "RR does not accept multiplier override for {}", mode));
        }
    }
}

auto BayesRRScheme::make_prior(const BayesModel& model) const
    -> std::vector<GeneticPrior>
{
    std::vector<GeneticPrior> priors;
    for (const auto mode : options_.modes)
    {
        const auto& effect
            = mode == GeneticMode::A ? options_.additive : options_.dominance;
        double h2 = 0.0;
        switch (mode)
        {
            case GeneticMode::A:
                h2 = effect.heritability() ? effect.heritability()->value()
                                           : 0.5;
                break;
            case GeneticMode::D:
                h2 = effect.heritability() ? effect.heritability()->value()
                                           : 0.2;
                break;
        }

        const auto* genetic = model.genetic(mode);
        if (genetic == nullptr)
        {
            throw GelexException(
                fmt::format("genetic design not found for mode {}", mode));
        }
        if (genetic->design_variance <= 0)
        {
            throw GelexException(
                fmt::format(
                    "genetic design_variance must be positive for mode {}",
                    mode));
        }

        const double active_marker_weight = 1.0;
        if (active_marker_weight <= 0)
        {
            throw GelexException(
                fmt::format(
                    "active_marker_weight must be positive, got {}",
                    active_marker_weight));
        }
        const double target
            = model.phenotype_variance() * h2
              / (active_marker_weight * genetic->design_variance);
        if (!std::isfinite(target) || target <= 0)
        {
            throw GelexException(
                fmt::format(
                    "target_marker_variance must be finite and positive, got "
                    "{}",
                    target));
        }

        priors.emplace_back(
            SingleGeneticPrior{SingleSharedGaussianPrior{
                mode,
                SharedMarkerVariance{VarianceParameter{
                    target,
                    ScaledInvChiSqPrior{4.0, (4.0 - 2.0) / 4.0 * target}}}}});
    }
    return priors;
}

BayesAScheme::BayesAScheme(const BayesRecipeConfig& options) : options_{options}
{
    if (options_.joint_proportion)
    {
        throw GelexException("A does not accept --joint-pi");
    }
    if (options_.joint_proportion_update)
    {
        throw GelexException("A does not accept --estimate-joint-pi");
    }
    for (const auto mode : options_.modes)
    {
        const auto& effect
            = mode == GeneticMode::A ? options_.additive : options_.dominance;
        if (effect.proportion() || effect.proportion_update())
        {
            throw GelexException(
                fmt::format(
                    "A does not accept proportion override for {}", mode));
        }
        if (effect.multiplier())
        {
            throw GelexException(
                fmt::format(
                    "A does not accept multiplier override for {}", mode));
        }
    }
}

auto BayesAScheme::make_prior(const BayesModel& model) const
    -> std::vector<GeneticPrior>
{
    std::vector<GeneticPrior> priors;
    for (const auto mode : options_.modes)
    {
        const auto& effect
            = mode == GeneticMode::A ? options_.additive : options_.dominance;
        double h2 = 0.0;
        switch (mode)
        {
            case GeneticMode::A:
                h2 = effect.heritability() ? effect.heritability()->value()
                                           : 0.5;
                break;
            case GeneticMode::D:
                h2 = effect.heritability() ? effect.heritability()->value()
                                           : 0.2;
                break;
        }

        const auto* genetic = model.genetic(mode);
        if (genetic == nullptr)
        {
            throw GelexException(
                fmt::format("genetic design not found for mode {}", mode));
        }
        if (genetic->design_variance <= 0)
        {
            throw GelexException(
                fmt::format(
                    "genetic design_variance must be positive for mode {}",
                    mode));
        }

        const double active_marker_weight = 1.0;
        if (active_marker_weight <= 0)
        {
            throw GelexException(
                fmt::format(
                    "active_marker_weight must be positive, got {}",
                    active_marker_weight));
        }
        const double target
            = model.phenotype_variance() * h2
              / (active_marker_weight * genetic->design_variance);
        if (!std::isfinite(target) || target <= 0)
        {
            throw GelexException(
                fmt::format(
                    "target_marker_variance must be finite and positive, got "
                    "{}",
                    target));
        }

        priors.emplace_back(
            SingleGeneticPrior{SinglePerMarkerGaussianPrior{
                mode,
                PerMarkerVariance{VarianceParameter{
                    target,
                    ScaledInvChiSqPrior{4.0, (4.0 - 2.0) / 4.0 * target}}}}});
    }
    return priors;
}

BayesBScheme::BayesBScheme(const BayesRecipeConfig& options) : options_{options}
{
    if (options_.joint_proportion)
    {
        throw GelexException("B does not accept --joint-pi");
    }
    if (options_.joint_proportion_update)
    {
        throw GelexException("B does not accept --estimate-joint-pi");
    }
    for (const auto mode : options_.modes)
    {
        const auto& effect
            = mode == GeneticMode::A ? options_.additive : options_.dominance;
        if (effect.multiplier())
        {
            throw GelexException(
                fmt::format(
                    "B does not accept multiplier override for {}", mode));
        }
    }
}

auto BayesBScheme::make_prior(const BayesModel& model) const
    -> std::vector<GeneticPrior>
{
    std::vector<GeneticPrior> priors;
    for (const auto mode : options_.modes)
    {
        const auto& effect
            = mode == GeneticMode::A ? options_.additive : options_.dominance;
        const Simplex<double> proportion
            = effect.proportion().value_or(Simplex<double>{{0.99, 0.01}});
        const double active_marker_weight = 1.0 - proportion[0];

        double h2 = 0.0;
        switch (mode)
        {
            case GeneticMode::A:
                h2 = effect.heritability() ? effect.heritability()->value()
                                           : 0.5;
                break;
            case GeneticMode::D:
                h2 = effect.heritability() ? effect.heritability()->value()
                                           : 0.2;
                break;
        }

        const auto* genetic = model.genetic(mode);
        if (genetic == nullptr)
        {
            throw GelexException(
                fmt::format("genetic design not found for mode {}", mode));
        }
        if (genetic->design_variance <= 0)
        {
            throw GelexException(
                fmt::format(
                    "genetic design_variance must be positive for mode {}",
                    mode));
        }
        if (active_marker_weight <= 0)
        {
            throw GelexException(
                fmt::format(
                    "active_marker_weight must be positive, got {}",
                    active_marker_weight));
        }
        const double target
            = model.phenotype_variance() * h2
              / (active_marker_weight * genetic->design_variance);
        if (!std::isfinite(target) || target <= 0)
        {
            throw GelexException(
                fmt::format(
                    "target_marker_variance must be finite and positive, got "
                    "{}",
                    target));
        }

        if (effect.proportion_update().value_or(true))
        {
            const auto n = static_cast<Eigen::Index>(proportion.size());
            priors.emplace_back(
                SingleGeneticPrior{SinglePerMarkerSpikeSlabGaussianPrior{
                    mode,
                    PerMarkerVariance{VarianceParameter{
                        target,
                        ScaledInvChiSqPrior{4.0, (4.0 - 2.0) / 4.0 * target}}},
                    MixtureProportion{SimplexParameter{
                        proportion.to_mat(),
                        DirichletPrior{Eigen::VectorXd::Ones(n)}}}}});
        }
        else
        {
            priors.emplace_back(
                SingleGeneticPrior{SinglePerMarkerSpikeSlabGaussianPrior{
                    mode,
                    PerMarkerVariance{VarianceParameter{
                        target,
                        ScaledInvChiSqPrior{4.0, (4.0 - 2.0) / 4.0 * target}}},
                    MixtureProportion{proportion.to_mat()}}});
        }
    }
    return priors;
}

BayesCScheme::BayesCScheme(const BayesRecipeConfig& options) : options_{options}
{
    if (options_.joint_proportion)
    {
        throw GelexException("C does not accept --joint-pi");
    }
    if (options_.joint_proportion_update)
    {
        throw GelexException("C does not accept --estimate-joint-pi");
    }
    for (const auto mode : options_.modes)
    {
        const auto& effect
            = mode == GeneticMode::A ? options_.additive : options_.dominance;
        if (effect.multiplier())
        {
            throw GelexException(
                fmt::format(
                    "C does not accept multiplier override for {}", mode));
        }
    }
}

auto BayesCScheme::make_prior(const BayesModel& model) const
    -> std::vector<GeneticPrior>
{
    std::vector<GeneticPrior> priors;
    for (const auto mode : options_.modes)
    {
        const auto& effect
            = mode == GeneticMode::A ? options_.additive : options_.dominance;
        const Simplex<double> proportion
            = effect.proportion().value_or(Simplex<double>{{0.99, 0.01}});
        const double active_marker_weight = 1.0 - proportion[0];

        double h2 = 0.0;
        switch (mode)
        {
            case GeneticMode::A:
                h2 = effect.heritability() ? effect.heritability()->value()
                                           : 0.5;
                break;
            case GeneticMode::D:
                h2 = effect.heritability() ? effect.heritability()->value()
                                           : 0.2;
                break;
        }

        const auto* genetic = model.genetic(mode);
        if (genetic == nullptr)
        {
            throw GelexException(
                fmt::format("genetic design not found for mode {}", mode));
        }
        if (genetic->design_variance <= 0)
        {
            throw GelexException(
                fmt::format(
                    "genetic design_variance must be positive for mode {}",
                    mode));
        }
        if (active_marker_weight <= 0)
        {
            throw GelexException(
                fmt::format(
                    "active_marker_weight must be positive, got {}",
                    active_marker_weight));
        }
        const double target
            = model.phenotype_variance() * h2
              / (active_marker_weight * genetic->design_variance);
        if (!std::isfinite(target) || target <= 0)
        {
            throw GelexException(
                fmt::format(
                    "target_marker_variance must be finite and positive, got "
                    "{}",
                    target));
        }

        if (effect.proportion_update().value_or(true))
        {
            const auto n = static_cast<Eigen::Index>(proportion.size());
            priors.emplace_back(
                SingleGeneticPrior{SingleSharedSpikeSlabGaussianPrior{
                    mode,
                    SharedMarkerVariance{VarianceParameter{
                        target,
                        ScaledInvChiSqPrior{4.0, (4.0 - 2.0) / 4.0 * target}}},
                    MixtureProportion{SimplexParameter{
                        proportion.to_mat(),
                        DirichletPrior{Eigen::VectorXd::Ones(n)}}}}});
        }
        else
        {
            priors.emplace_back(
                SingleGeneticPrior{SingleSharedSpikeSlabGaussianPrior{
                    mode,
                    SharedMarkerVariance{VarianceParameter{
                        target,
                        ScaledInvChiSqPrior{4.0, (4.0 - 2.0) / 4.0 * target}}},
                    MixtureProportion{proportion.to_mat()}}});
        }
    }
    return priors;
}

BayesRScheme::BayesRScheme(const BayesRecipeConfig& options) : options_{options}
{
    if (options_.joint_proportion)
    {
        throw GelexException("R does not accept --joint-pi");
    }
    if (options_.joint_proportion_update)
    {
        throw GelexException("R does not accept --estimate-joint-pi");
    }
    for (const auto mode : options_.modes)
    {
        const auto& effect
            = mode == GeneticMode::A ? options_.additive : options_.dominance;
        if (effect.proportion().has_value() != effect.multiplier().has_value())
        {
            throw GelexException(
                fmt::format(
                    "R: proportion and multiplier must be paired for {}",
                    mode));
        }
    }
}

auto BayesRScheme::make_prior(const BayesModel& model) const
    -> std::vector<GeneticPrior>
{
    std::vector<GeneticPrior> priors;
    for (const auto mode : options_.modes)
    {
        const auto& effect
            = mode == GeneticMode::A ? options_.additive : options_.dominance;
        const Simplex<double> proportion = effect.proportion().value_or(
            Simplex<double>{{0.99, 0.005, 0.003, 0.001, 0.001}});
        const ScaleMultiplier<double> multiplier = effect.multiplier().value_or(
            ScaleMultiplier<double>{{0.0, 0.001, 0.01, 0.1, 1.0}});
        double active_marker_weight = 0.0;
        for (std::size_t i = 0; i < proportion.size(); ++i)
        {
            active_marker_weight += proportion[i] * multiplier[i];
        }

        double h2 = 0.0;
        switch (mode)
        {
            case GeneticMode::A:
                h2 = effect.heritability() ? effect.heritability()->value()
                                           : 0.5;
                break;
            case GeneticMode::D:
                h2 = effect.heritability() ? effect.heritability()->value()
                                           : 0.2;
                break;
        }

        const auto* genetic = model.genetic(mode);
        if (genetic == nullptr)
        {
            throw GelexException(
                fmt::format("genetic design not found for mode {}", mode));
        }
        if (genetic->design_variance <= 0)
        {
            throw GelexException(
                fmt::format(
                    "genetic design_variance must be positive for mode {}",
                    mode));
        }
        if (active_marker_weight <= 0)
        {
            throw GelexException(
                fmt::format(
                    "active_marker_weight must be positive, got {}",
                    active_marker_weight));
        }
        const double target
            = model.phenotype_variance() * h2
              / (active_marker_weight * genetic->design_variance);
        if (!std::isfinite(target) || target <= 0)
        {
            throw GelexException(
                fmt::format(
                    "target_marker_variance must be finite and positive, got "
                    "{}",
                    target));
        }

        if (effect.proportion_update().value_or(true))
        {
            const auto n = static_cast<Eigen::Index>(proportion.size());
            priors.emplace_back(
                SingleGeneticPrior{SingleScaledMixtureGaussianPrior{
                    mode,
                    SharedMarkerVariance{VarianceParameter{
                        target,
                        ScaledInvChiSqPrior{4.0, (4.0 - 2.0) / 4.0 * target}}},
                    multiplier.to_mat(),
                    MixtureProportion{SimplexParameter{
                        proportion.to_mat(),
                        DirichletPrior{Eigen::VectorXd::Ones(n)}}}}});
        }
        else
        {
            priors.emplace_back(
                SingleGeneticPrior{SingleScaledMixtureGaussianPrior{
                    mode,
                    SharedMarkerVariance{VarianceParameter{
                        target,
                        ScaledInvChiSqPrior{4.0, (4.0 - 2.0) / 4.0 * target}}},
                    multiplier.to_mat(),
                    MixtureProportion{proportion.to_mat()}}});
        }
    }
    return priors;
}

}  // namespace gelex::bayes
