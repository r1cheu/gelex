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
#include <cmath>
#include <cstddef>
#include <random>
#include <type_traits>

#include "gelex/bayes/stats/scaled_inv_chi2_distribution.h"

using Catch::Approx;
using gelex::ScaledInvChi2Distribution;

TEST_CASE(
    "ScaledInvChi2Distribution follows the parameter API",
    "[bayes][stats][scaled_inv_chi2_distribution]")
{
    using Distribution = ScaledInvChi2Distribution<>;
    using Parameters = Distribution::param_type;
    static_assert(std::is_same_v<Parameters::distribution_type, Distribution>);

    const Parameters default_parameters;
    REQUIRE(default_parameters.degrees_of_freedom() == 1.0);
    REQUIRE(default_parameters.scale() == 1.0);

    const Parameters parameters{4.0, 2.0};
    Distribution distribution;
    distribution.param(parameters);
    REQUIRE(distribution.param() == parameters);
    REQUIRE(distribution.degrees_of_freedom() == 4.0);
    REQUIRE(distribution.scale() == 2.0);
    REQUIRE(Distribution{parameters}.param() == parameters);
    REQUIRE(Distribution{4.0, 2.0}.param() == parameters);
    REQUIRE(distribution.min() == 0.0);
    REQUIRE(std::isinf(distribution.max()));

    distribution.reset();
    REQUIRE(distribution.param() == parameters);

    std::mt19937_64 stored_rng{481};
    std::mt19937_64 supplied_rng{481};
    Distribution supplied_distribution;
    for (std::size_t draw_index = 0; draw_index < 100; ++draw_index)
    {
        REQUIRE(
            distribution(stored_rng)
            == supplied_distribution(supplied_rng, parameters));
    }
}

TEST_CASE(
    "ScaledInvChi2Distribution samples the requested parameters",
    "[bayes][stats][scaled_inv_chi2_distribution]")
{
    constexpr std::size_t draw_count = 50'000;
    constexpr double degrees_of_freedom = 30.0;
    constexpr double scale = 2.0;
    ScaledInvChi2Distribution distribution{degrees_of_freedom, scale};
    std::mt19937_64 rng{123};
    double mean = 0.0;
    double sum_squared_deviations = 0.0;
    bool all_positive = true;

    for (std::size_t draw_index = 0; draw_index < draw_count; ++draw_index)
    {
        const double sample = distribution(rng);
        all_positive = all_positive && sample > 0.0;
        const double delta = sample - mean;
        mean += delta / static_cast<double>(draw_index + 1);
        sum_squared_deviations += delta * (sample - mean);
    }

    const double variance
        = sum_squared_deviations / static_cast<double>(draw_count - 1);
    const double expected_mean
        = (degrees_of_freedom * scale) / (degrees_of_freedom - 2.0);
    const double expected_variance
        = (2.0 * degrees_of_freedom * degrees_of_freedom * scale * scale)
          / ((degrees_of_freedom - 2.0) * (degrees_of_freedom - 2.0)
             * (degrees_of_freedom - 4.0));
    REQUIRE(all_positive);
    REQUIRE(mean == Approx(expected_mean).margin(0.02));
    REQUIRE(variance == Approx(expected_variance).margin(0.02));
}
