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
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "gelex/exception.h"
#include "gelex/infra/stats/detail/normal.h"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using gelex::detail::norm_cdf;
using gelex::detail::norm_ppf;

// Reference values from scipy.stats.norm.ppf
TEST_CASE("norm_ppf known values", "[math_utils]")
{
    SECTION("median")
    {
        CHECK_THAT(norm_ppf(0.5), WithinAbs(0.0, 1e-15));
    }

    SECTION("standard quantiles")
    {
        // scipy.stats.norm.ppf reference values
        CHECK_THAT(
            norm_ppf(0.025), WithinRel(-1.959963984540054494e+00, 1e-12));
        CHECK_THAT(norm_ppf(0.975), WithinRel(1.959963984540054049e+00, 1e-12));
        CHECK_THAT(norm_ppf(0.05), WithinRel(-1.644853626951472858e+00, 1e-12));
        CHECK_THAT(norm_ppf(0.95), WithinRel(1.644853626951472192e+00, 1e-12));
        CHECK_THAT(norm_ppf(0.1), WithinRel(-1.281551565544600368e+00, 1e-12));
        CHECK_THAT(norm_ppf(0.9), WithinRel(1.281551565544600368e+00, 1e-12));
    }

    SECTION("sigma boundaries")
    {
        CHECK_THAT(norm_ppf(8.413447460685429258e-01), WithinRel(1.0, 1e-12));
        CHECK_THAT(norm_ppf(1.586552539314570742e-01), WithinRel(-1.0, 1e-12));
        CHECK_THAT(norm_ppf(9.772498680518207914e-01), WithinRel(2.0, 1e-12));
        CHECK_THAT(norm_ppf(9.986501019683698965e-01), WithinRel(3.0, 1e-12));
    }

    SECTION("extreme tails")
    {
        CHECK_THAT(
            norm_ppf(0.001), WithinRel(-3.090232306167813192e+00, 1e-12));
        CHECK_THAT(norm_ppf(0.999), WithinRel(3.090232306167813192e+00, 1e-12));
        CHECK_THAT(norm_ppf(1e-6), WithinRel(-4.753424308822898681e+00, 1e-10));
        CHECK_THAT(
            norm_ppf(1.0 - 1e-6), WithinRel(4.753424308822898681e+00, 1e-10));
        CHECK_THAT(
            norm_ppf(1e-10), WithinRel(-6.361340902404055697e+00, 1e-10));
    }
}

TEST_CASE("norm_ppf boundary throws", "[math_utils]")
{
    CHECK_THROWS_AS(norm_ppf(0.0), gelex::GelexException);
    CHECK_THROWS_AS(norm_ppf(1.0), gelex::GelexException);
    CHECK_THROWS_AS(norm_ppf(-0.1), gelex::GelexException);
    CHECK_THROWS_AS(norm_ppf(1.1), gelex::GelexException);
}

TEST_CASE("norm_ppf symmetry", "[math_utils]")
{
    for (double p : {0.01, 0.05, 0.1, 0.2, 0.3, 0.4})
    {
        CHECK_THAT(norm_ppf(p), WithinAbs(-norm_ppf(1.0 - p), 1e-14));
    }
}

TEST_CASE("norm_ppf roundtrip with norm_cdf", "[math_utils]")
{
    for (double p :
         {0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999})
    {
        double x = norm_ppf(p);
        double roundtrip = norm_cdf(x);
        CHECK_THAT(roundtrip, WithinAbs(p, 1e-14));
    }
}

TEST_CASE("norm_ppf accuracy vs norm_cdf", "[math_utils]")
{
    // Sweep across probabilities and verify norm_cdf(norm_ppf(p)) ≈ p
    for (int i = 1; i < 1000; ++i)
    {
        double p = static_cast<double>(i) / 1000.0;
        double x = norm_ppf(p);
        double recovered = norm_cdf(x);
        CHECK_THAT(recovered, WithinAbs(p, 1e-13));
    }
}
