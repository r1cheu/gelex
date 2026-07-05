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

#include "bayes/detail/scheme_helpers.h"

#include <cmath>
#include <cstddef>
#include <optional>
#include <string_view>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/recipe_options.h"
#include "gelex/exception.h"
#include "gelex/types/constrained_vector.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::bayes::detail
{

auto heritability(const BayesRecipeOptions& options, GeneticMode mode) -> double
{
    switch (mode)
    {
        case GeneticMode::A:
            return options.additive_heritability
                       ? options.additive_heritability->value()
                       : 0.5;
        case GeneticMode::D:
            return options.dominance_heritability
                       ? options.dominance_heritability->value()
                       : 0.2;
    }
    throw GelexException(fmt::format("Unsupported genetic mode: {}", mode));
}

auto proportion(const BayesRecipeOptions& options, GeneticMode mode)
    -> const std::optional<Simplex<double>>&
{
    switch (mode)
    {
        case GeneticMode::A:
            return options.additive_proportion;
        case GeneticMode::D:
            return options.dominance_proportion;
    }
    throw GelexException(fmt::format("Unsupported genetic mode: {}", mode));
}

auto multiplier(const BayesRecipeOptions& options, GeneticMode mode)
    -> const std::optional<ScaleMultiplier<double>>&
{
    switch (mode)
    {
        case GeneticMode::A:
            return options.additive_multiplier;
        case GeneticMode::D:
            return options.dominance_multiplier;
    }
    throw GelexException(fmt::format("Unsupported genetic mode: {}", mode));
}

auto proportion_update(const BayesRecipeOptions& options, GeneticMode mode)
    -> const std::optional<bool>&
{
    switch (mode)
    {
        case GeneticMode::A:
            return options.additive_proportion_update;
        case GeneticMode::D:
            return options.dominance_proportion_update;
    }
    throw GelexException(fmt::format("Unsupported genetic mode: {}", mode));
}

auto target_marker_variance(
    const BayesModel& model,
    GeneticMode mode,
    double h2,
    double active_marker_weight) -> double
{
    const auto* genetic = model.genetic(mode);
    if (genetic == nullptr)
    {
        throw GelexException(
            fmt::format("genetic design not found for mode {}", mode));
    }
    const double design_variance = genetic->col_var.sum();
    if (design_variance <= 0)
    {
        throw GelexException(
            fmt::format(
                "genetic col_var sum must be positive for mode {}", mode));
    }
    if (active_marker_weight <= 0)
    {
        throw GelexException(
            fmt::format(
                "active_marker_weight must be positive, got {}",
                active_marker_weight));
    }

    const double target = model.phenotype_variance() * h2
                          / (active_marker_weight * design_variance);
    if (!std::isfinite(target) || target <= 0)
    {
        throw GelexException(
            fmt::format(
                "target_marker_variance must be finite and positive, got {}",
                target));
    }
    return target;
}

auto variance_parameter(double target) -> VarianceParameter
{
    return VarianceParameter{
        target, ScaledInvChiSqPrior{4.0, (4.0 - 2.0) / 4.0 * target}};
}

auto mixture_proportion(const Simplex<double>& proportion, bool sampled)
    -> MixtureProportion
{
    if (sampled)
    {
        const auto n = static_cast<Eigen::Index>(proportion.size());
        return MixtureProportion{SimplexParameter{
            proportion.to_mat(), DirichletPrior{Eigen::VectorXd::Ones(n)}}};
    }
    return MixtureProportion{proportion.to_mat()};
}

auto scaled_active_marker_weight(
    const Simplex<double>& proportion,
    const ScaleMultiplier<double>& multiplier) -> double
{
    double active_marker_weight = 0.0;
    for (std::size_t i = 0; i < proportion.size(); ++i)
    {
        active_marker_weight += proportion[i] * multiplier[i];
    }
    return active_marker_weight;
}

auto reject_joint_options(
    const BayesRecipeOptions& options,
    std::string_view scheme) -> void
{
    if (options.joint_proportion)
    {
        throw GelexException(fmt::format("{} does not accept --jpi", scheme));
    }
    if (options.joint_proportion_update)
    {
        throw GelexException(
            fmt::format("{} does not accept --sample-jpi", scheme));
    }
    if (options.dominance_positive_probability)
    {
        throw GelexException(
            fmt::format("{} does not accept --dom-pos-prob", scheme));
    }
}

auto reject_proportion_options(
    const BayesRecipeOptions& options,
    std::string_view scheme) -> void
{
    for (const auto mode : options.modes)
    {
        if (proportion(options, mode) || proportion_update(options, mode))
        {
            throw GelexException(
                fmt::format(
                    "{} does not accept proportion override for {}",
                    scheme,
                    mode));
        }
    }
}

auto reject_multiplier_options(
    const BayesRecipeOptions& options,
    std::string_view scheme) -> void
{
    for (const auto mode : options.modes)
    {
        if (multiplier(options, mode))
        {
            throw GelexException(
                fmt::format(
                    "{} does not accept multiplier override for {}",
                    scheme,
                    mode));
        }
    }
}

auto reject_unpaired_proportion_multiplier(
    const BayesRecipeOptions& options,
    std::string_view scheme) -> void
{
    for (const auto mode : options.modes)
    {
        const auto& mode_proportion = proportion(options, mode);
        const auto& mode_multiplier = multiplier(options, mode);
        if (mode_proportion.has_value() != mode_multiplier.has_value())
        {
            throw GelexException(
                fmt::format(
                    "{}: proportion and multiplier must be paired for {}",
                    scheme,
                    mode));
        }
        if (mode_proportion && mode_multiplier
            && mode_proportion->size() != mode_multiplier->size())
        {
            throw GelexException(
                fmt::format(
                    "{}: proportion size {} != multiplier size {} for {}",
                    scheme,
                    mode_proportion->size(),
                    mode_multiplier->size(),
                    mode));
        }
    }
}

}  // namespace gelex::bayes::detail
