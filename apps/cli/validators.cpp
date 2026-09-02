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

#include "cli/validators.h"

#include <CLI/CLI.hpp>
#include <string>
#include <vector>

#include "gelex/bayes/builtin_method.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/rank_inverse_norm_transform.h"
#include "gelex/genetic_mode.h"

namespace cli
{

auto open_unit_interval() -> CLI::Validator
{
    return CLI::Validator{
        [](std::string& input)
        {
            double value{};
            if (!CLI::detail::lexical_cast(input, value) || value <= 0.0
                || value >= 1.0)
            {
                return std::string{
                    "expected a value in the open interval (0, 1)"};
            }
            return std::string{};
        },
        "FLOAT in (0 - 1)",
        "PROB"};
}

auto non_negative_number() -> CLI::Validator
{
    return CLI::Validator{
        [](std::string& input)
        {
            double value{};
            if (!CLI::detail::lexical_cast(input, value) || !(value >= 0.0))
            {
                return std::string{"expected a non-negative value"};
            }
            return std::string{};
        },
        "NONNEGATIVE",
        "NONNEGATIVE"};
}

auto genotype_method_validator() -> CLI::Validator
{
    std::vector<std::string> names;
    names.reserve(gelex::GENOTYPE_METHOD_CODES.size());
    for (const auto& entry : gelex::GENOTYPE_METHOD_CODES)
    {
        names.emplace_back(entry.first);
    }
    return CLI::IsMember(names, CLI::ignore_case);
}

auto bayes_method_validator() -> CLI::Validator
{
    std::vector<std::string> names;
    names.reserve(gelex::bayes_method_names.size());
    for (const auto& [method, name] : gelex::bayes_method_names)
    {
        names.emplace_back(name);
    }
    return CLI::IsMember(names);
}

auto genetic_mode_set_validator() -> CLI::Validator
{
    std::vector<std::string> names;
    names.reserve(gelex::genetic_mode_set_names.size());
    for (const auto& [set, name] : gelex::genetic_mode_set_names)
    {
        names.emplace_back(name);
    }
    return CLI::IsMember(names);
}

auto rint_type_validator() -> CLI::Validator
{
    std::vector<std::string> names;
    names.reserve(gelex::rint_type_names.size());
    for (const auto& [type, name] : gelex::rint_type_names)
    {
        names.emplace_back(name);
    }
    return CLI::IsMember(names);
}

}  // namespace cli
