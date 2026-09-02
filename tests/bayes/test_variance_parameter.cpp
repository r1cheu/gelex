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
#include <limits>

#include "gelex/bayes/variance_parameter.h"
#include "gelex/exception.h"

TEST_CASE(
    "scaled inverse chi-squared prior accepts proper parameters",
    "[bayes][variance_parameter]")
{
    const auto prior = gelex::ScaledInvChiSqPrior{4.0, 1.5};

    REQUIRE(prior.degrees_of_freedom() == 4.0);
    REQUIRE(prior.scale() == 1.5);
}

TEST_CASE(
    "scaled inverse chi-squared prior rejects improper parameters",
    "[bayes][variance_parameter]")
{
    SECTION("non-positive degrees of freedom")
    {
        REQUIRE_THROWS_AS(
            (gelex::ScaledInvChiSqPrior{0.0, 1.0}), gelex::GelexException);
    }

    SECTION("non-positive scale")
    {
        REQUIRE_THROWS_AS(
            (gelex::ScaledInvChiSqPrior{4.0, 0.0}), gelex::GelexException);
    }

    SECTION("non-finite values")
    {
        const double infinity = std::numeric_limits<double>::infinity();
        REQUIRE_THROWS_AS(
            (gelex::ScaledInvChiSqPrior{infinity, 1.0}), gelex::GelexException);
        REQUIRE_THROWS_AS(
            (gelex::ScaledInvChiSqPrior{4.0, infinity}), gelex::GelexException);
    }
}

TEST_CASE(
    "variance parameter requires a positive finite initial value",
    "[bayes][variance_parameter]")
{
    const auto prior = gelex::ScaledInvChiSqPrior{4.0, 1.0};

    REQUIRE_THROWS_AS(
        (gelex::VarianceParameter{0.0, prior}), gelex::GelexException);
    REQUIRE_THROWS_AS(
        (gelex::VarianceParameter{
            std::numeric_limits<double>::infinity(), prior}),
        gelex::GelexException);
}
