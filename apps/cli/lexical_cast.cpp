// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "cli/lexical_cast.h"

#include <algorithm>
#include <cctype>
#include <string>

#include "gelex/bayes/builtin_method.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/rank_inverse_norm_transform.h"
#include "gelex/genetic_mode.h"

namespace gelex
{

auto lexical_cast(const std::string& input, BayesMethod& output) -> bool
{
    for (const auto& [method, name] : bayes_method_names)
    {
        if (input == name)
        {
            output = method;
            return true;
        }
    }
    return false;
}

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
