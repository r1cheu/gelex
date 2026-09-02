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

#include "config.h"

#include <fmt/format.h>
#include <optional>
#include <string_view>
#include <vector>

#include "gelex/bayes/builtin_method.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

namespace cli
{

namespace
{

struct ModeOptionNames
{
    std::string_view additive;
    std::string_view dominance;
};

template <typename T>
auto is_supplied(const std::optional<T>& value) noexcept -> bool
{
    return value.has_value();
}

template <typename T, typename Allocator>
auto is_supplied(const std::vector<T, Allocator>& values) noexcept -> bool
{
    return !values.empty();
}

template <gelex::GeneticMode Mode>
constexpr auto option_name(ModeOptionNames names) noexcept -> std::string_view
{
    if constexpr (Mode == gelex::GeneticMode::A)
    {
        return names.additive;
    }
    else
    {
        static_assert(Mode == gelex::GeneticMode::D);
        return names.dominance;
    }
}

auto method_name(gelex::BayesMethod method) noexcept -> std::string_view
{
    for (const auto& [value, name] : gelex::bayes_method_names)
    {
        if (value == method)
        {
            return name;
        }
    }
    return "unknown";
}

template <typename Values>
auto validate_mode_options(
    const Values& values,
    ModeOptionNames names,
    const McmcConfig& config,
    bool method_accepts_options) -> void
{
    values.for_each(
        [&]<gelex::GeneticMode Mode>(const auto& value)
        {
            if (!is_supplied(value))
            {
                return;
            }

            const auto option = option_name<Mode>(names);
            if (!config.mode.contains(Mode))
            {
                throw gelex::GelexException(
                    fmt::format(
                        "{} requires --mode to include {}", option, Mode));
            }
            if (!method_accepts_options)
            {
                throw gelex::GelexException(
                    fmt::format(
                        "{} is not valid for --method {}",
                        option,
                        method_name(config.method)));
            }
        });
}

constexpr auto accepts_mode_probabilities(gelex::BayesMethod method) noexcept
    -> bool
{
    return method == gelex::BayesMethod::B || method == gelex::BayesMethod::C
           || method == gelex::BayesMethod::R;
}

constexpr auto accepts_mode_scales(gelex::BayesMethod method) noexcept -> bool
{
    return method == gelex::BayesMethod::R;
}

}  // namespace

auto validate_mcmc_config(const McmcConfig& config) -> void
{
    if (config.method == gelex::BayesMethod::CD
        && config.mode != McmcConfig::option_modes)
    {
        throw gelex::GelexException("--method CD requires --mode AD");
    }

    validate_mode_options(
        config.genetic_variance_shares,
        {.additive = "--h2", .dominance = "--d2"},
        config,
        true);
    validate_mode_options(
        config.mixture_probabilities.mode_values(),
        {.additive = "--pi", .dominance = "--dpi"},
        config,
        accepts_mode_probabilities(config.method));
    validate_mode_options(
        config.mixture_scales,
        {.additive = "--scale", .dominance = "--dscale"},
        config,
        accepts_mode_scales(config.method));

    if (!config.mixture_probabilities.joint().empty()
        && config.method != gelex::BayesMethod::CD)
    {
        throw gelex::GelexException(
            fmt::format(
                "--jpi is not valid for --method {}",
                method_name(config.method)));
    }

    if (config.dominance_positive_probability)
    {
        if (!config.mode.contains(gelex::GeneticMode::D))
        {
            throw gelex::GelexException(
                "--dom-pos-prob requires --mode to include D");
        }
        if (config.method != gelex::BayesMethod::CD)
        {
            throw gelex::GelexException(
                fmt::format(
                    "--dom-pos-prob is not valid for --method {}",
                    method_name(config.method)));
        }
    }

    if (config.random.has_random_design() != config.random_pve.has_value())
    {
        throw gelex::GelexException(
            "--random-pve and at least one of --drand/--qrand must be "
            "provided together");
    }
}

}  // namespace cli
