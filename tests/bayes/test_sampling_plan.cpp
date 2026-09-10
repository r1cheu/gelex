// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <catch2/catch_test_macros.hpp>

#include "gelex/bayes/sampling_plan.h"
#include "gelex/exception.h"

TEST_CASE("Sampling plan retains post burn-in draws", "[bayes][sampling_plan]")
{
    const gelex::SamplingPlan plan{5, 1, 2, 42};
    REQUIRE(plan.draw_count() == 2);
    REQUIRE(plan.seed() == 42);
    REQUIRE_FALSE(plan.retains(0));
    REQUIRE_FALSE(plan.retains(1));
    REQUIRE(plan.retains(2));
    REQUIRE_FALSE(plan.retains(3));
    REQUIRE(plan.retains(4));
}

TEST_CASE("Sampling plan rejects invalid schedules", "[bayes][sampling_plan]")
{
    REQUIRE_THROWS_AS(gelex::SamplingPlan(0, 0, 1, 0), gelex::GelexException);
    REQUIRE_THROWS_AS(gelex::SamplingPlan(-1, 0, 1, 0), gelex::GelexException);
    REQUIRE_THROWS_AS(gelex::SamplingPlan(4, -1, 1, 0), gelex::GelexException);
    REQUIRE_THROWS_AS(gelex::SamplingPlan(4, 4, 1, 0), gelex::GelexException);
    REQUIRE_THROWS_AS(gelex::SamplingPlan(4, 0, 0, 0), gelex::GelexException);
    REQUIRE_THROWS_AS(gelex::SamplingPlan(4, 0, -1, 0), gelex::GelexException);
    REQUIRE_THROWS_AS(gelex::SamplingPlan(4, 1, 2, 0), gelex::GelexException);
}
