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
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <concepts>
#include <limits>
#include <string>

#include "gelex/bayes/variance/budget.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using gelex::default_additive_share;
using gelex::default_shares;
using gelex::GelexException;
using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::VarianceBudget;

namespace
{

// Shares are a single input object: they cannot be supplied one at a time, and
// there is no empty budget standing for "unspecified".
static_assert(!std::default_initializable<VarianceBudget>);
static_assert(!std::constructible_from<VarianceBudget, double, double, double>);

// residual() is a sum and only ever agrees with a decimal literal up to
// rounding, so it is checked at run time with a matcher; share() returns what
// was stored and is compared exactly.
constexpr double residual_tolerance = 1e-12;

const auto budget
    = VarianceBudget{{.additive = 0.4, .dominance = 0.05, .random = 0.05}};

const auto additive_only
    = VarianceBudget{default_shares(GeneticModeSet{GeneticMode::A})};
const auto dominance_only
    = VarianceBudget{default_shares(GeneticModeSet{GeneticMode::D})};
const auto both
    = VarianceBudget{default_shares(GeneticMode::A | GeneticMode::D)};

// An input adapter overrides the fields the user gave and keeps the rest, so
// the untouched share holds its default without the caller naming it. Asserted
// against the constant rather than a literal: this pins the override mechanism,
// not the value of the table.
const auto partially_overridden = []
{
    auto shares = default_shares(GeneticMode::A | GeneticMode::D);
    shares.dominance = 0.3;
    return VarianceBudget{shares};
}();

auto message_of(auto&& construct) -> std::string
{
    try
    {
        construct();
    }
    catch (const GelexException& error)
    {
        return error.what();
    }
    return {};
}

}  // namespace

TEST_CASE(
    "VarianceBudget derives the residual from every share",
    "[bayes][variance_budget]")
{
    REQUIRE(budget.share(GeneticMode::A) == 0.4);
    REQUIRE(budget.share(GeneticMode::D) == 0.05);
    REQUIRE(budget.random() == 0.05);
    REQUIRE_THAT(budget.residual(), WithinAbs(0.5, residual_tolerance));
    REQUIRE_THAT(additive_only.residual(), WithinAbs(0.5, residual_tolerance));
    REQUIRE_THAT(dominance_only.residual(), WithinAbs(0.8, residual_tolerance));
    REQUIRE_THAT(both.residual(), WithinAbs(0.3, residual_tolerance));
}

TEST_CASE("VarianceBudget rejects invalid shares", "[bayes][variance_budget]")
{
    SECTION("a share is negative")
    {
        REQUIRE_THAT(
            message_of([] { return VarianceBudget{{.dominance = -0.1}}; }),
            ContainsSubstring(
                "dominance variance share must be finite and non-negative"));
    }

    SECTION("a share is not finite")
    {
        REQUIRE_THAT(
            message_of(
                []
                {
                    return VarianceBudget{
                        {.random = std::numeric_limits<double>::quiet_NaN()}};
                }),
            ContainsSubstring(
                "random variance share must be finite and non-negative"));
    }

    SECTION("positive infinity is not finite")
    {
        REQUIRE_THAT(
            message_of(
                []
                {
                    return VarianceBudget{
                        {.additive = std::numeric_limits<double>::infinity()}};
                }),
            ContainsSubstring(
                "additive variance share must be finite and non-negative"));
    }

    SECTION("negative infinity is not finite")
    {
        REQUIRE_THAT(
            message_of(
                []
                {
                    return VarianceBudget{
                        {.random = -std::numeric_limits<double>::infinity()}};
                }),
            ContainsSubstring(
                "random variance share must be finite and non-negative"));
    }

    SECTION("shares consume the residual")
    {
        REQUIRE_THAT(
            message_of(
                []
                {
                    return VarianceBudget{{.additive = 0.6, .dominance = 0.4}};
                }),
            ContainsSubstring("variance shares must sum to less than 1"));
    }
}

TEST_CASE(
    "VarianceBudget accepts every strictly positive residual share",
    "[bayes][variance_budget]")
{
    REQUIRE_NOTHROW(
        VarianceBudget{
            {.additive = std::numeric_limits<double>::denorm_min()}});

    const auto budget_at_upper_boundary
        = VarianceBudget{{.additive = std::nextafter(1.0, 0.0)}};
    REQUIRE(budget_at_upper_boundary.residual() > 0.0);
}

TEST_CASE(
    "default_shares gives absent modes a zero share",
    "[bayes][variance_budget]")
{
    REQUIRE(
        partially_overridden.share(GeneticMode::A) == default_additive_share);
    REQUIRE(partially_overridden.share(GeneticMode::D) == 0.3);

    for (const auto modes : {
             GeneticModeSet{GeneticMode::A},
             GeneticModeSet{GeneticMode::D},
             GeneticMode::A | GeneticMode::D,
         })
    {
        const auto budget = VarianceBudget{default_shares(modes)};

        for (const auto mode : gelex::all_genetic_modes)
        {
            REQUIRE((budget.share(mode) > 0.0) == modes.contains(mode));
        }
    }
}
