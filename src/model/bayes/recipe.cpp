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

#include "gelex/model/bayes/recipe.h"

#include <memory>
#include <string_view>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "independent_recipes.h"
#include "joint_recipes.h"
#include "recipe_impl.h"

namespace gelex::bayes
{

BayesRecipe::BayesRecipe(BayesRecipePreset preset, BayesRecipeConfig options)
    : preset_{preset}, options_{std::move(options)}
{
    validate_modes(options_.modes);
    impl_ = make_impl(preset_, options_);
}

BayesRecipe::~BayesRecipe() = default;

auto BayesRecipe::make_prior(const BayesModel& model) const -> BayesPrior
{
    return BayesPrior{
        make_random_prior(model),
        impl_->make_genetic_prior_blocks(model),
        make_residual_prior(model)};
}

auto BayesRecipe::make_impl(
    BayesRecipePreset preset,
    const BayesRecipeConfig& options) -> std::unique_ptr<BayesRecipeImpl>
{
    switch (preset)
    {
        case BayesRecipePreset::RR:
            return std::make_unique<BayesRRMethod>(options);
        case BayesRecipePreset::A:
            return std::make_unique<BayesAMethod>(options);
        case BayesRecipePreset::B:
            return std::make_unique<BayesBMethod>(options);
        case BayesRecipePreset::C:
            return std::make_unique<BayesCMethod>(options);
        case BayesRecipePreset::R:
            return std::make_unique<BayesRMethod>(options);
        case BayesRecipePreset::CD:
            return std::make_unique<BayesCDMethod>(options);
        default:
            throw GelexException(
                fmt::format("Unsupported BayesRecipePreset: {}", preset));
    }
}

auto BayesRecipe::validate_modes(std::span<const GeneticMode> modes) -> void
{
    const bool valid
        = (modes.size() == 1
           && (modes[0] == GeneticMode::A || modes[0] == GeneticMode::D))
          || (modes.size() == 2 && modes[0] == GeneticMode::A
              && modes[1] == GeneticMode::D);
    if (!valid)
    {
        throw GelexException(
            "BayesRecipeConfig modes must be {A}, {D}, or {A, D}");
    }
}

auto BayesRecipe::make_random_prior(const BayesModel& model) const
    -> RandomPrior
{
    const double proportion = options_.random_variance_proportion
                                  ? options_.random_variance_proportion->value()
                                  : 0.05;
    return RandomPrior(
        model.phenotype_variance() * proportion, ScaledInvChiSqPrior{-2, 0});
}

auto BayesRecipe::make_residual_prior(const BayesModel& model) -> ResidualPrior
{
    const double default_residual_proportion = 0.5;
    return ResidualPrior(
        model.phenotype_variance() * default_residual_proportion,
        ScaledInvChiSqPrior{-2, 0});
}

auto to_bayes_recipe_preset(std::string_view preset) -> BayesRecipePreset
{
    for (const auto& [value, name] : kBayesRecipePresetNames)
    {
        if (preset == name)
        {
            return value;
        }
    }
    throw GelexException(fmt::format("Unknown method preset: {}", preset));
}

}  // namespace gelex::bayes
