// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <limits>

#include "gelex/bayes/variance/budget.h"
#include "gelex/genetic_mode.h"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using gelex::default_additive_share;
using gelex::default_dominance_share;
using gelex::default_proportion;
using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::VarianceBudget;

TEST_CASE(
    "VarianceBudget derives the residual from allocated proportions",
    "[bayes][variance_budget]")
{
    constexpr double tolerance = 1e-12;
    const auto budget
        = VarianceBudget{{.additive = 0.4, .dominance = 0.05, .random = 0.05}};

    REQUIRE(budget.genetic(GeneticMode::A) == 0.4);
    REQUIRE(budget.genetic(GeneticMode::D) == 0.05);
    REQUIRE(budget.random() == 0.05);
    REQUIRE_THAT(budget.residual(), WithinAbs(0.5, tolerance));
}

TEST_CASE("VarianceBudget rejects invalid shares", "[bayes][variance_budget]")
{
    SECTION("a share is negative")
    {
        REQUIRE_THROWS_WITH(
            (VarianceBudget{{.dominance = -0.1}}),
            ContainsSubstring(
                "dominance variance share must be finite and non-negative"));
    }

    SECTION("a share is not finite")
    {
        REQUIRE_THROWS_WITH(
            (VarianceBudget{
                {.random = std::numeric_limits<double>::quiet_NaN()}}),
            ContainsSubstring(
                "random variance share must be finite and non-negative"));
    }

    SECTION("positive infinity is not finite")
    {
        REQUIRE_THROWS_WITH(
            (VarianceBudget{
                {.additive = std::numeric_limits<double>::infinity()}}),
            ContainsSubstring(
                "additive variance share must be finite and non-negative"));
    }

    SECTION("shares consume the residual")
    {
        REQUIRE_THROWS_WITH(
            (VarianceBudget{{.additive = 0.6, .dominance = 0.4}}),
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
    "default_proportion gives absent modes a zero proportion",
    "[bayes][variance_budget]")
{
    const auto additive_only
        = default_proportion(GeneticModeSet{GeneticMode::A});
    REQUIRE(additive_only.additive == default_additive_share);
    REQUIRE(additive_only.dominance == 0.0);

    const auto dominance_only
        = default_proportion(GeneticModeSet{GeneticMode::D});
    REQUIRE(dominance_only.additive == 0.0);
    REQUIRE(dominance_only.dominance == default_dominance_share);

    auto proportions = default_proportion(GeneticMode::A | GeneticMode::D);
    proportions.dominance = 0.3;
    const auto partially_overridden = VarianceBudget{proportions};

    REQUIRE(
        partially_overridden.genetic(GeneticMode::A) == default_additive_share);
    REQUIRE(partially_overridden.genetic(GeneticMode::D) == 0.3);
}
