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

#include <algorithm>
#include <array>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <CLI/CLI.hpp>

#include "cli/cli_helper.h"
#include "gelex/bayes/recipe.h"
#include "gelex/exception.h"
#include "gelex/types/constrained_value.h"
#include "gelex/types/constrained_vector.h"
#include "gelex/types/genetic_effect_type.h"

namespace cli
{

namespace
{

auto has_mode(
    std::span<const gelex::GeneticMode> modes,
    gelex::GeneticMode mode) -> bool
{
    return std::ranges::any_of(
        modes,
        [mode](gelex::GeneticMode candidate) { return candidate == mode; });
}

template <typename T, typename Raw = T>
auto get_optional(const CLI::App& cmd, std::string_view arg) -> std::optional<T>
{
    auto* option = cmd.get_option(std::string{arg});
    if (option->count() == 0)
    {
        return std::nullopt;
    }
    return T{option->as<Raw>()};
}

auto reject_effect_flags_without_mode(
    const CLI::App& cmd,
    std::span<const gelex::GeneticMode> modes,
    gelex::GeneticMode mode,
    std::span<const std::string_view> flags) -> void
{
    if (has_mode(modes, mode))
    {
        return;
    }
    for (const auto flag : flags)
    {
        if (cmd.get_option(std::string{flag})->count() > 0)
        {
            throw gelex::GelexException(
                fmt::format("{} requires --mode to include {}", flag, mode));
        }
    }
}

auto reject_effect_flags_without_mode(
    const CLI::App& cmd,
    std::span<const gelex::GeneticMode> modes) -> void
{
    constexpr std::array ADDITIVE_FLAGS
        = {std::string_view{"--h2"},
           std::string_view{"--pi"},
           std::string_view{"--scale"},
           std::string_view{"--sample-pi"}};
    constexpr std::array DOMINANCE_FLAGS
        = {std::string_view{"--d2"},
           std::string_view{"--dom-pos-prob"},
           std::string_view{"--dpi"},
           std::string_view{"--dscale"},
           std::string_view{"--sample-dpi"}};

    reject_effect_flags_without_mode(
        cmd, modes, gelex::GeneticMode::A, ADDITIVE_FLAGS);
    reject_effect_flags_without_mode(
        cmd, modes, gelex::GeneticMode::D, DOMINANCE_FLAGS);
}

}  // namespace

auto make_bayes_recipe_options(const CLI::App& cmd)
    -> gelex::bayes::BayesRecipeOptions
{
    auto modes
        = parse_genetic_modes(cmd.get_option("--mode")->as<std::string>());
    reject_effect_flags_without_mode(cmd, modes);

    return gelex::bayes::BayesRecipeOptions{
        .scheme = gelex::bayes::to_bayes_recipe_scheme(
            cmd.get_option("--method")->as<std::string>()),
        .modes = std::move(modes),
        .additive_heritability
        = get_optional<gelex::OpenUnitInterval<double>, double>(cmd, "--h2"),
        .additive_proportion
        = get_optional<gelex::Simplex<double>, std::vector<double>>(
            cmd, "--pi"),
        .additive_multiplier
        = get_optional<gelex::ScaleMultiplier<double>, std::vector<double>>(
            cmd, "--scale"),
        .additive_proportion_update = get_optional<bool>(cmd, "--sample-pi"),
        .dominance_heritability
        = get_optional<gelex::OpenUnitInterval<double>, double>(cmd, "--d2"),
        .dominance_proportion
        = get_optional<gelex::Simplex<double>, std::vector<double>>(
            cmd, "--dpi"),
        .dominance_multiplier
        = get_optional<gelex::ScaleMultiplier<double>, std::vector<double>>(
            cmd, "--dscale"),
        .dominance_proportion_update = get_optional<bool>(cmd, "--sample-dpi"),
        .dominance_positive_probability
        = get_optional<gelex::OpenUnitInterval<double>, double>(
            cmd, "--dom-pos-prob"),
        .joint_proportion
        = get_optional<gelex::Simplex<double>, std::vector<double>>(
            cmd, "--jpi"),
        .joint_proportion_update = get_optional<bool>(cmd, "--sample-jpi"),
        .random_variance_proportion
        = get_optional<gelex::OpenUnitInterval<double>, double>(
            cmd, "--random-pve"),
    };
}

}  // namespace cli
