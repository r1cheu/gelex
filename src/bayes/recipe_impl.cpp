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

#include "recipe_impl.h"

#include <cmath>
#include <string_view>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe_options.h"
#include "gelex/exception.h"
#include "gelex/types/constrained_value.h"
#include "gelex/types/constrained_vector.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

BayesRecipeImpl::BayesRecipeImpl(
    std::string_view name,
    const BayesRecipeConfig& options)
    : name_{name}, options_{options}
{
}

auto BayesRecipeImpl::default_heritability(GeneticMode mode)
    -> OpenUnitInterval<double>
{
    switch (mode)
    {
        case GeneticMode::A:
            return OpenUnitInterval<double>{0.5};
        case GeneticMode::D:
            return OpenUnitInterval<double>{0.2};
    }
    std::unreachable();
}

auto BayesRecipeImpl::marker_variance_from_heritability(
    const BayesModel& model,
    GeneticMode mode,
    double heritability,
    double active_marker_weight) -> double
{
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
                "genetic design_variance must be positive for mode {}", mode));
    }
    if (active_marker_weight <= 0)
    {
        throw GelexException(
            fmt::format(
                "active_marker_weight must be positive, got {}",
                active_marker_weight));
    }

    const double target = model.phenotype_variance() * heritability
                          / (active_marker_weight * genetic->design_variance);

    if (!std::isfinite(target) || target <= 0)
    {
        throw GelexException(
            fmt::format(
                "target_marker_variance must be finite and positive, got {}",
                target));
    }

    return target;
}

auto BayesRecipeImpl::reject_dominance_positive_probability_override() const
    -> void
{
    if (options_.dominance_positive_probability)
    {
        throw GelexException(
            fmt::format(
                "{} does not support --dominance-positive-prob", name_));
    }
}

IndependentMethod::IndependentMethod(
    std::string_view name,
    const BayesRecipeConfig& options)
    : BayesRecipeImpl{name, options}
{
}

auto IndependentMethod::make_genetic_prior_blocks(const BayesModel& model) const
    -> std::vector<GeneticPrior>
{
    std::vector<GeneticPrior> priors;
    for_each_effect(
        [&](GeneticMode mode, const EffectConfig& effect)
        {
            priors.emplace_back(make_single_genetic_prior(mode, effect, model));
        });
    return priors;
}

auto IndependentMethod::reject_joint_overrides() const -> void
{
    if (options().joint_proportion)
    {
        throw GelexException(
            fmt::format("{} does not accept --joint-pi", name()));
    }
    if (options().joint_proportion_update)
    {
        throw GelexException(
            fmt::format("{} does not accept --estimate-joint-pi", name()));
    }
}

auto IndependentMethod::reject_per_effect_proportion() const -> void
{
    for_each_effect(
        [&](GeneticMode mode, const EffectConfig& effect)
        {
            if (effect.proportion() || effect.proportion_update())
            {
                throw GelexException(
                    fmt::format(
                        "{} does not accept proportion override for {}",
                        name(),
                        mode));
            }
        });
}

auto IndependentMethod::reject_per_effect_multiplier() const -> void
{
    for_each_effect(
        [&](GeneticMode mode, const EffectConfig& effect)
        {
            if (effect.multiplier())
            {
                throw GelexException(
                    fmt::format(
                        "{} does not accept multiplier override for {}",
                        name(),
                        mode));
            }
        });
}

auto IndependentMethod::require_paired_proportion_and_multiplier() const -> void
{
    for_each_effect(
        [&](GeneticMode mode, const EffectConfig& effect)
        {
            if (effect.proportion().has_value()
                != effect.multiplier().has_value())
            {
                throw GelexException(
                    fmt::format(
                        "{}: proportion and multiplier must be paired for {}",
                        name(),
                        mode));
            }
        });
}

JointMethod::JointMethod(
    std::string_view name,
    const BayesRecipeConfig& options)
    : BayesRecipeImpl{name, options}
{
}

auto JointMethod::make_genetic_prior_blocks(const BayesModel& model) const
    -> std::vector<GeneticPrior>
{
    std::vector<GeneticPrior> priors;
    priors.emplace_back(make_joint_prior(options(), model));
    return priors;
}

auto JointMethod::require_both_modes() const -> void
{
    if (options().modes.size() != 2)
    {
        throw GelexException(
            fmt::format("{} requires --mode AD (both effects)", name()));
    }
}

auto JointMethod::reject_per_effect_proportion() const -> void
{
    auto check = [&](GeneticMode mode, const EffectConfig& effect)
    {
        if (effect.proportion() || effect.proportion_update())
        {
            throw GelexException(
                fmt::format(
                    "{} does not accept per-effect proportion override for {}; "
                    "use --joint-pi instead",
                    name(),
                    mode));
        }
    };
    check(GeneticMode::A, options().additive);
    check(GeneticMode::D, options().dominance);
}

auto JointMethod::reject_per_effect_multiplier() const -> void
{
    auto check = [&](GeneticMode mode, const EffectConfig& effect)
    {
        if (effect.multiplier())
        {
            throw GelexException(
                fmt::format(
                    "{} does not accept per-effect multiplier override for {}",
                    name(),
                    mode));
        }
    };
    check(GeneticMode::A, options().additive);
    check(GeneticMode::D, options().dominance);
}

}  // namespace gelex::bayes
