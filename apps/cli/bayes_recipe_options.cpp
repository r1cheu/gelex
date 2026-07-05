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

#include "bayes_recipe_options.h"

#include <optional>
#include <string_view>

#include <fmt/format.h>

#include "gelex/bayes/recipe.h"
#include "gelex/exception.h"
#include "gelex/types/constrained_value.h"
#include "gelex/types/constrained_vector.h"
#include "gelex/types/genetic_mode.h"

namespace cli
{

namespace
{

auto reject_effect_flags_without_mode(
    const McmcConfig& config,
    gelex::GeneticModeSet modes) -> void
{
    const auto require_mode
        = [&](bool present, std::string_view flag, gelex::GeneticMode mode)
    {
        if (present && !modes.contains(mode))
        {
            throw gelex::GelexException(
                fmt::format("{} requires --mode to include {}", flag, mode));
        }
    };
    require_mode(config.h2.has_value(), "--h2", gelex::GeneticMode::A);
    require_mode(!config.pi.empty(), "--pi", gelex::GeneticMode::A);
    require_mode(!config.scale.empty(), "--scale", gelex::GeneticMode::A);
    require_mode(
        config.sample_pi.has_value(), "--sample-pi", gelex::GeneticMode::A);
    require_mode(config.d2.has_value(), "--d2", gelex::GeneticMode::D);
    require_mode(
        config.dom_pos_prob.has_value(),
        "--dom-pos-prob",
        gelex::GeneticMode::D);
    require_mode(!config.dpi.empty(), "--dpi", gelex::GeneticMode::D);
    require_mode(!config.dscale.empty(), "--dscale", gelex::GeneticMode::D);
    require_mode(
        config.sample_dpi.has_value(), "--sample-dpi", gelex::GeneticMode::D);
}

}  // namespace

auto make_bayes_recipe_options(const McmcConfig& config)
    -> gelex::bayes::BayesRecipeOptions
{
    reject_effect_flags_without_mode(config, config.mode);

    return gelex::bayes::BayesRecipeOptions{
        .scheme = gelex::bayes::to_bayes_recipe_scheme(config.method),
        .modes = config.mode,
        .additive_heritability
        = config.h2 ? std::optional{gelex::OpenUnitInterval<double>{*config.h2}}
                    : std::nullopt,
        .additive_proportion
        = config.pi.empty() ? std::nullopt
                            : std::optional{gelex::Simplex<double>{config.pi}},
        .additive_multiplier
        = config.scale.empty()
              ? std::nullopt
              : std::optional{gelex::ScaleMultiplier<double>{config.scale}},
        .additive_proportion_update = config.sample_pi,
        .dominance_heritability
        = config.d2 ? std::optional{gelex::OpenUnitInterval<double>{*config.d2}}
                    : std::nullopt,
        .dominance_proportion
        = config.dpi.empty()
              ? std::nullopt
              : std::optional{gelex::Simplex<double>{config.dpi}},
        .dominance_multiplier
        = config.dscale.empty()
              ? std::nullopt
              : std::optional{gelex::ScaleMultiplier<double>{config.dscale}},
        .dominance_proportion_update = config.sample_dpi,
        .dominance_positive_probability
        = config.dom_pos_prob ? std::optional{gelex::OpenUnitInterval<double>{
                                    *config.dom_pos_prob}}
                              : std::nullopt,
        .joint_proportion
        = config.jpi.empty()
              ? std::nullopt
              : std::optional{gelex::Simplex<double>{config.jpi}},
        .joint_proportion_update = config.sample_jpi,
        .random_variance_proportion
        = config.random_pve ? std::optional{gelex::OpenUnitInterval<double>{
                                  *config.random_pve}}
                            : std::nullopt,
    };
}

}  // namespace cli
