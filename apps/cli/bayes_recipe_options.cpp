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

#include <argparse.h>
#include <fmt/format.h>

#include "cli/cli_helper.h"
#include "gelex/bayes/recipe.h"
#include "gelex/exception.h"
#include "gelex/types/constrained_value.h"
#include "gelex/types/constrained_vector.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::cli
{

namespace
{

auto has_mode(std::span<const GeneticMode> modes, GeneticMode mode) -> bool
{
    return std::ranges::any_of(
        modes, [mode](GeneticMode candidate) { return candidate == mode; });
}

template <typename T, typename Raw = T>
auto get_optional(const argparse::ArgumentParser& cmd, std::string_view arg)
    -> std::optional<T>
{
    if (!cmd.is_used(arg))
    {
        return std::nullopt;
    }
    return T{cmd.get<Raw>(arg)};
}

auto reject_effect_flags_without_mode(
    const argparse::ArgumentParser& cmd,
    std::span<const GeneticMode> modes,
    GeneticMode mode,
    std::span<const std::string_view> flags) -> void
{
    if (has_mode(modes, mode))
    {
        return;
    }
    for (const auto flag : flags)
    {
        if (cmd.is_used(flag))
        {
            throw GelexException(
                fmt::format("{} requires --mode to include {}", flag, mode));
        }
    }
}

auto reject_effect_flags_without_mode(
    const argparse::ArgumentParser& cmd,
    std::span<const GeneticMode> modes) -> void
{
    constexpr std::array ADDITIVE_FLAGS
        = {std::string_view{"--additive-h2"},
           std::string_view{"--additive-pi"},
           std::string_view{"--additive-multiplier"},
           std::string_view{"--estimate-additive-pi"}};
    constexpr std::array DOMINANCE_FLAGS
        = {std::string_view{"--dominance-h2"},
           std::string_view{"--dominance-pi"},
           std::string_view{"--dominance-multiplier"},
           std::string_view{"--estimate-dominance-pi"}};

    reject_effect_flags_without_mode(
        cmd, modes, GeneticMode::A, ADDITIVE_FLAGS);
    reject_effect_flags_without_mode(
        cmd, modes, GeneticMode::D, DOMINANCE_FLAGS);
}

}  // namespace

auto make_bayes_recipe_options(const argparse::ArgumentParser& cmd)
    -> bayes::BayesRecipeOptions
{
    auto modes = parse_genetic_modes(cmd.get<std::string>("--mode"));
    reject_effect_flags_without_mode(cmd, modes);

    return bayes::BayesRecipeOptions{
        .scheme
        = bayes::to_bayes_recipe_scheme(cmd.get<std::string>("--method")),
        .modes = std::move(modes),
        .additive_heritability
        = get_optional<gelex::OpenUnitInterval<double>, double>(
            cmd, "--additive-h2"),
        .additive_proportion
        = get_optional<gelex::Simplex<double>, std::vector<double>>(
            cmd, "--additive-pi"),
        .additive_multiplier
        = get_optional<gelex::ScaleMultiplier<double>, std::vector<double>>(
            cmd, "--additive-multiplier"),
        .additive_proportion_update
        = get_optional<bool>(cmd, "--estimate-additive-pi"),
        .dominance_heritability
        = get_optional<gelex::OpenUnitInterval<double>, double>(
            cmd, "--dominance-h2"),
        .dominance_proportion
        = get_optional<gelex::Simplex<double>, std::vector<double>>(
            cmd, "--dominance-pi"),
        .dominance_multiplier
        = get_optional<gelex::ScaleMultiplier<double>, std::vector<double>>(
            cmd, "--dominance-multiplier"),
        .dominance_proportion_update
        = get_optional<bool>(cmd, "--estimate-dominance-pi"),
        .joint_proportion
        = get_optional<gelex::Simplex<double>, std::vector<double>>(
            cmd, "--joint-pi"),
        .joint_proportion_update
        = get_optional<bool>(cmd, "--estimate-joint-pi"),
        .random_variance_proportion
        = get_optional<gelex::OpenUnitInterval<double>, double>(
            cmd, "--random-variance-proportion"),
    };
}

}  // namespace gelex::cli
