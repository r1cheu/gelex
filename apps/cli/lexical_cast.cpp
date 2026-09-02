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

#include "cli/lexical_cast.h"

#include <algorithm>
#include <cctype>
#include <string>

#include "gelex/data/genotype_method.h"
#include "gelex/data/rank_inverse_norm_transform.h"
#include "gelex/genetic_mode.h"

namespace gelex
{

auto lexical_cast(const std::string& input, GenotypeMethod& output) -> bool
{
    for (const auto& [code, method] : GENOTYPE_METHOD_CODES)
    {
        if (input.size() == code.size()
            && std::equal(
                code.begin(),
                code.end(),
                input.begin(),
                [](unsigned char expected, unsigned char actual)
                { return std::tolower(expected) == std::tolower(actual); }))
        {
            output = method;
            return true;
        }
    }
    return false;
}

auto lexical_cast(const std::string& input, GeneticMode& output) -> bool
{
    for (const auto& [mode, name] : genetic_mode_names)
    {
        if (input == name)
        {
            output = mode;
            return true;
        }
    }
    return false;
}

auto lexical_cast(const std::string& input, GeneticModeSet& output) -> bool
{
    for (const auto& [set, name] : genetic_mode_set_names)
    {
        if (input == name)
        {
            output = set;
            return true;
        }
    }
    return false;
}

auto lexical_cast(const std::string& input, RintType& output) -> bool
{
    for (const auto& [type, name] : rint_type_names)
    {
        if (input == name)
        {
            output = type;
            return true;
        }
    }
    return false;
}

}  // namespace gelex
