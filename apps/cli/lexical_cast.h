// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_LEXICAL_CAST_H_
#define APPS_CLI_LEXICAL_CAST_H_

#include <CLI/Error.hpp>
#include <string>

#include "gelex/bayes/builtin_method.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/rank_inverse_norm_transform.h"
#include "gelex/genetic_mode.h"

// Found by CLI11 via ADL, so add_option can bind enum targets directly.
namespace gelex
{
auto lexical_cast(const std::string& input, BayesMethod& output) -> bool;

auto lexical_cast(const std::string& input, GenotypeMethod& output) -> bool;

auto lexical_cast(const std::string& input, GeneticMode& output) -> bool;

auto lexical_cast(const std::string& input, GeneticModeSet& output) -> bool;

auto lexical_cast(const std::string& input, RintType& output) -> bool;
}  // namespace gelex

namespace cli
{

// add_option's generic path assigns a default-constructed target for empty
// input. Types whose values are all meaningful have no such default, so they
// bind through add_option_function and this callback instead.
template <typename T>
auto lexical_assigner(T& target)
{
    return [&target](const std::string& input)
    {
        if (!gelex::lexical_cast(input, target))
        {
            throw CLI::ValidationError("cannot parse '" + input + "'");
        }
    };
}

}  // namespace cli

#endif  // APPS_CLI_LEXICAL_CAST_H_
