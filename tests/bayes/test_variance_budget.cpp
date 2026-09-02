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
#include <concepts>

#include "gelex/bayes/variance_budget.h"
#include "gelex/types/genetic_mode.h"

using Catch::Matchers::WithinAbs;
using gelex::DEFAULT_ADDITIVE_SHARE;
using gelex::default_shares;
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
constexpr double RESIDUAL_TOLERANCE = 1e-12;

constexpr auto BUDGET
    = VarianceBudget{{.additive = 0.4, .dominance = 0.05, .random = 0.05}};

static_assert(BUDGET.share(GeneticMode::A) == 0.4);
static_assert(BUDGET.share(GeneticMode::D) == 0.05);
static_assert(BUDGET.random() == 0.05);

constexpr auto ADDITIVE_ONLY
    = VarianceBudget{default_shares(GeneticModeSet{GeneticMode::A})};
constexpr auto DOMINANCE_ONLY
    = VarianceBudget{default_shares(GeneticModeSet{GeneticMode::D})};
constexpr auto BOTH
    = VarianceBudget{default_shares(GeneticMode::A | GeneticMode::D)};

static_assert(ADDITIVE_ONLY.share(GeneticMode::A) == 0.5);
static_assert(ADDITIVE_ONLY.share(GeneticMode::D) == 0.0);

static_assert(DOMINANCE_ONLY.share(GeneticMode::A) == 0.0);
static_assert(DOMINANCE_ONLY.share(GeneticMode::D) == 0.2);

static_assert(BOTH.share(GeneticMode::A) == 0.5);
static_assert(BOTH.share(GeneticMode::D) == 0.2);

// Defaults never request a random design.
static_assert(ADDITIVE_ONLY.random() == 0.0);
static_assert(DOMINANCE_ONLY.random() == 0.0);
static_assert(BOTH.random() == 0.0);

// An input adapter overrides the fields the user gave and keeps the rest, so
// the untouched share holds its default without the caller naming it. Asserted
// against the constant rather than a literal: this pins the override mechanism,
// not the value of the table.
constexpr auto PARTIALLY_OVERRIDDEN = []
{
    auto shares = default_shares(GeneticMode::A | GeneticMode::D);
    shares.dominance = 0.3;
    return VarianceBudget{shares};
}();

static_assert(
    PARTIALLY_OVERRIDDEN.share(GeneticMode::A) == DEFAULT_ADDITIVE_SHARE);
static_assert(PARTIALLY_OVERRIDDEN.share(GeneticMode::D) == 0.3);

}  // namespace

TEST_CASE(
    "VarianceBudget derives the residual from every share",
    "[bayes][variance_budget]")
{
    REQUIRE_THAT(BUDGET.residual(), WithinAbs(0.5, RESIDUAL_TOLERANCE));
    REQUIRE_THAT(ADDITIVE_ONLY.residual(), WithinAbs(0.5, RESIDUAL_TOLERANCE));
    REQUIRE_THAT(DOMINANCE_ONLY.residual(), WithinAbs(0.8, RESIDUAL_TOLERANCE));
    REQUIRE_THAT(BOTH.residual(), WithinAbs(0.3, RESIDUAL_TOLERANCE));
}

TEST_CASE(
    "VarianceBudget represents an over-allocated budget",
    "[bayes][variance_budget]")
{
    // Construction does not validate, so this is representable and reports a
    // negative residual. SpecDiagnostics rejects it later.
    REQUIRE_THAT(
        VarianceBudget{{.additive = 1.5}}.residual(),
        WithinAbs(-0.5, RESIDUAL_TOLERANCE));
}

TEST_CASE(
    "default_shares gives absent modes a zero share",
    "[bayes][variance_budget]")
{
    for (const auto modes : {
             GeneticModeSet{GeneticMode::A},
             GeneticModeSet{GeneticMode::D},
             GeneticMode::A | GeneticMode::D,
         })
    {
        const auto budget = VarianceBudget{default_shares(modes)};

        for (const auto mode : gelex::ALL_GENETIC_MODES)
        {
            REQUIRE((budget.share(mode) > 0.0) == modes.contains(mode));
        }
    }
}
