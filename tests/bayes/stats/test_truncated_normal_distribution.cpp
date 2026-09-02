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
#include <cstdint>
#include <limits>
#include <numbers>
#include <random>

#include "gelex/bayes/stats/truncated_normal_distribution.h"

namespace gelex
{
namespace
{

using Catch::Approx;

constexpr std::uint64_t seed = 0xDEADBEEF42ULL;

TEST_CASE(
    "TruncatedNormalDistribution follows the STL parameter API",
    "[bayes][stats][truncated_normal_distribution]")
{
    using Distribution = TruncatedNormalDistribution<>;

    const Distribution::param_type default_parameters;
    REQUIRE(default_parameters.mean() == 0.0);
    REQUIRE(default_parameters.stddev() == 1.0);
    REQUIRE(default_parameters.support() == HalfLine::Positive);

    const Distribution::param_type parameters{-1.5, 2.0, HalfLine::Negative};
    Distribution distribution{parameters};
    REQUIRE(distribution.mean() == -1.5);
    REQUIRE(distribution.stddev() == 2.0);
    REQUIRE(distribution.support() == HalfLine::Negative);
    REQUIRE(distribution.param() == parameters);
    REQUIRE(distribution.min() == std::numeric_limits<double>::lowest());
    REQUIRE(distribution.max() == 0.0);

    const Distribution::param_type replacement{0.5, 0.75};
    distribution.param(replacement);
    REQUIRE(distribution.param() == replacement);
    REQUIRE(distribution.min() == 0.0);
    REQUIRE(distribution.max() == std::numeric_limits<double>::max());
}

TEST_CASE(
    "TruncatedNormalDistribution samples the selected half line",
    "[bayes][stats][truncated_normal_distribution]")
{
    TruncatedNormalDistribution distribution;
    std::mt19937_64 rng{seed};

    const TruncatedNormalDistribution<>::param_type positive{
        1.0, 2.0, HalfLine::Positive};
    const TruncatedNormalDistribution<>::param_type negative{
        -1.0, 2.0, HalfLine::Negative};
    for (std::size_t draw_index = 0; draw_index < 10'000; ++draw_index)
    {
        REQUIRE(distribution(rng, positive) > 0.0);
        REQUIRE(distribution(rng, negative) < 0.0);
    }
}

TEST_CASE(
    "TruncatedNormalDistribution matches standard half-normal moments",
    "[bayes][stats][truncated_normal_distribution]")
{
    constexpr std::size_t draw_count = 200'000;
    TruncatedNormalDistribution distribution;
    std::mt19937_64 rng{seed};

    double sum = 0.0;
    double sum_squares = 0.0;
    for (std::size_t draw_index = 0; draw_index < draw_count; ++draw_index)
    {
        const double value = distribution(rng);
        sum += value;
        sum_squares += value * value;
    }

    const double mean = sum / static_cast<double>(draw_count);
    const double variance
        = (sum_squares / static_cast<double>(draw_count)) - (mean * mean);
    REQUIRE(mean == Approx(std::sqrt(2.0 / std::numbers::pi)).margin(0.01));
    REQUIRE(variance == Approx(1.0 - (2.0 / std::numbers::pi)).margin(0.01));
}

TEST_CASE(
    "TruncatedNormalDistribution samples a remote tail",
    "[bayes][stats][truncated_normal_distribution]")
{
    constexpr std::size_t draw_count = 200'000;
    constexpr double mean = -3.0;
    constexpr double stddev = 1.0;
    constexpr double standardized_boundary = 3.0;
    TruncatedNormalDistribution distribution{mean, stddev};
    std::mt19937_64 rng{seed};

    double sum = 0.0;
    for (std::size_t draw_index = 0; draw_index < draw_count; ++draw_index)
    {
        sum += distribution(rng);
    }

    const double density
        = std::exp(-0.5 * standardized_boundary * standardized_boundary)
          / std::sqrt(2.0 * std::numbers::pi);
    const double tail_probability
        = 0.5 * std::erfc(standardized_boundary / std::numbers::sqrt2);
    const double expected_mean = mean + (stddev * density / tail_probability);
    REQUIRE(
        sum / static_cast<double>(draw_count)
        == Approx(expected_mean).margin(0.01));
}

TEST_CASE(
    "TruncatedNormalDistribution reset clears cached state",
    "[bayes][stats][truncated_normal_distribution]")
{
    TruncatedNormalDistribution distribution;
    std::mt19937_64 rng{seed};
    (void)distribution(rng);

    distribution.reset();
    rng.seed(seed);
    const double after_reset = distribution(rng);

    TruncatedNormalDistribution fresh_distribution;
    std::mt19937_64 fresh_rng{seed};
    REQUIRE(after_reset == fresh_distribution(fresh_rng));
}

}  // namespace
}  // namespace gelex
