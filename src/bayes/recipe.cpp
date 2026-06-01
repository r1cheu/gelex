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

#include "gelex/bayes/recipe.h"

#include <span>
#include <string_view>
#include <utility>
#include <variant>

#include <fmt/format.h>

#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe_options.h"
#include "gelex/exception.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

BayesRecipe::BayesRecipe(BayesRecipeOptions options)
    : options_{std::move(options)},
      scheme_{
          [this]() -> BayesScheme
          {
              validate_modes(options_.modes);
              switch (options_.scheme)
              {
                  case BayesRecipeScheme::RR:
                      return BayesRRScheme{options_};
                  case BayesRecipeScheme::A:
                      return BayesAScheme{options_};
                  case BayesRecipeScheme::B:
                      return BayesBScheme{options_};
                  case BayesRecipeScheme::C:
                      return BayesCScheme{options_};
                  case BayesRecipeScheme::R:
                      return BayesRScheme{options_};
                  case BayesRecipeScheme::CD:
                      return BayesCDScheme{options_};
                  default:
                      throw GelexException(
                          fmt::format(
                              "Unsupported BayesRecipeScheme: {}",
                              options_.scheme));
              }
          }()}
{
}

BayesRecipe::~BayesRecipe() = default;

auto BayesRecipe::make_prior(const BayesModel& model) const -> BayesPrior
{
    return BayesPrior{
        make_random_prior(model),
        std::visit(
            [&model](const auto& recipe) { return recipe.make_prior(model); },
            scheme_),
        make_residual_prior(model)};
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
            "BayesRecipeOptions modes must be {A}, {D}, or {A, D}");
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

auto to_bayes_recipe_scheme(std::string_view recipe_scheme) -> BayesRecipeScheme
{
    for (const auto& [value, name] : BAYES_RECIPE_SCHEME_NAMES)
    {
        if (recipe_scheme == name)
        {
            return value;
        }
    }
    throw GelexException(
        fmt::format("Unknown recipe scheme: {}", recipe_scheme));
}

}  // namespace gelex::bayes
