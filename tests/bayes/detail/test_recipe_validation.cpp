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

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <string>

#include "gelex/bayes/detail/recipe_validation.h"
#include "gelex/bayes/variance_budget.h"
#include "gelex/exception.h"
#include "gelex/types/genetic_mode.h"

using Catch::Matchers::ContainsSubstring;
using gelex::GelexException;
using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::VarianceBudget;
using gelex::detail::validate_recipe_inputs;

namespace
{

constexpr auto mode_a = GeneticModeSet{GeneticMode::A};
constexpr auto mode_ad = GeneticMode::A | GeneticMode::D;

const auto budget_a = VarianceBudget{{.additive = 0.5}};
const auto budget_ad = VarianceBudget{{.additive = 0.5, .dominance = 0.2}};

auto message_of(auto&& validate) -> std::string
{
    try
    {
        validate();
    }
    catch (const GelexException& error)
    {
        return error.what();
    }
    return {};
}

}  // namespace

TEST_CASE(
    "validate_recipe_inputs stays silent on valid input",
    "[bayes][recipe_validation]")
{
    REQUIRE_NOTHROW(validate_recipe_inputs(budget_a, mode_a));
    REQUIRE_NOTHROW(validate_recipe_inputs(budget_ad, mode_ad));
}

TEST_CASE(
    "validate_recipe_inputs cross-checks the budget against modes",
    "[bayes][recipe_validation]")
{
    SECTION("a present mode needs a positive share")
    {
        const auto message = message_of(
            []
            {
                validate_recipe_inputs(
                    VarianceBudget{{.additive = 0.5}}, mode_ad);
            });

        REQUIRE_THAT(
            message,
            ContainsSubstring(
                "D variance share must be positive when the mode is present"));
    }

    SECTION("an absent mode needs a zero share")
    {
        const auto message
            = message_of([] { validate_recipe_inputs(budget_ad, mode_a); });

        REQUIRE_THAT(
            message, ContainsSubstring("D variance share must be zero"));
    }
}
