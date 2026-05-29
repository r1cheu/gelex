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

#include "bayes_recipe_config.h"

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
#include "gelex/exception.h"
#include "gelex/model/bayes/recipe_options.h"
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

auto make_effect_config(
    const argparse::ArgumentParser& cmd,
    std::string_view h2_flag,
    std::string_view pi_flag,
    std::string_view multiplier_flag,
    std::string_view estimate_pi_flag) -> bayes::EffectConfig
{
    return bayes::EffectConfig{
        get_optional<gelex::OpenUnitInterval<double>, double>(cmd, h2_flag),
        get_optional<gelex::Simplex<double>, std::vector<double>>(cmd, pi_flag),
        get_optional<gelex::ScaleMultiplier<double>, std::vector<double>>(
            cmd, multiplier_flag),
        get_optional<bool>(cmd, estimate_pi_flag)};
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
    constexpr std::array kAdditiveFlags
        = {std::string_view{"--additive-h2"},
           std::string_view{"--additive-pi"},
           std::string_view{"--additive-multiplier"},
           std::string_view{"--estimate-additive-pi"}};
    constexpr std::array kDominanceFlags
        = {std::string_view{"--dominance-h2"},
           std::string_view{"--dominance-pi"},
           std::string_view{"--dominance-multiplier"},
           std::string_view{"--estimate-dominance-pi"},
           std::string_view{"--dominance-positive-prob"}};

    reject_effect_flags_without_mode(
        cmd, modes, GeneticMode::A, kAdditiveFlags);
    reject_effect_flags_without_mode(
        cmd, modes, GeneticMode::D, kDominanceFlags);
}

}  // namespace

auto make_bayes_recipe_config(const argparse::ArgumentParser& cmd)
    -> bayes::BayesRecipeConfig
{
    auto modes = parse_genetic_modes(cmd.get<std::string>("--mode"));
    reject_effect_flags_without_mode(cmd, modes);

    return bayes::BayesRecipeConfig{
        std::move(modes),
        make_effect_config(
            cmd,
            "--additive-h2",
            "--additive-pi",
            "--additive-multiplier",
            "--estimate-additive-pi"),
        make_effect_config(
            cmd,
            "--dominance-h2",
            "--dominance-pi",
            "--dominance-multiplier",
            "--estimate-dominance-pi"),
        get_optional<gelex::Simplex<double>, std::vector<double>>(
            cmd, "--joint-pi"),
        get_optional<bool>(cmd, "--estimate-joint-pi"),
        get_optional<gelex::OpenUnitInterval<double>, double>(
            cmd, "--dominance-positive-prob"),
        get_optional<gelex::OpenUnitInterval<double>, double>(
            cmd, "--random-variance-proportion")};
}

}  // namespace gelex::cli
