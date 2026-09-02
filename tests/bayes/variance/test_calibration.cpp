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

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <limits>

#include "gelex/bayes/variance/detail/calibration.h"
#include "gelex/exception.h"

using Catch::Approx;

TEST_CASE(
    "mean-calibrated variance parameter owns its initial value and prior",
    "[bayes][variance][calibration]")
{
    constexpr double target = 2.0;
    const auto parameter
        = gelex::detail::make_mean_calibrated_variance_parameter(target);
    const auto prior = parameter.prior.scaled_inv_chi2_parameters();

    REQUIRE(parameter.initial == target);
    REQUIRE(prior.degrees_of_freedom() == 4.0);
    REQUIRE(prior.scale() == Approx(0.5 * target));
}

TEST_CASE(
    "mean-calibrated variance parameter rejects invalid targets",
    "[bayes][variance][calibration]")
{
    SECTION("non-positive")
    {
        REQUIRE_THROWS_AS(
            gelex::detail::make_mean_calibrated_variance_parameter(0.0),
            gelex::GelexException);
    }

    SECTION("non-finite")
    {
        REQUIRE_THROWS_AS(
            gelex::detail::make_mean_calibrated_variance_parameter(
                std::numeric_limits<double>::infinity()),
            gelex::GelexException);
    }
}
