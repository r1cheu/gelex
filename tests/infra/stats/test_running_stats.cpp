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

#include <Eigen/Core>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <concepts>
#include <utility>

#include "gelex/exception.h"
#include "gelex/infra/stats/detail/running_stats.h"
#include "gelex/infra/stats/result.h"

namespace gelex
{

using Catch::Approx;
using gelex::detail::CategoricalFrequency;
using gelex::detail::ScalarRunningStats;
using gelex::detail::VectorRunningStats;

TEST_CASE(
    "ScalarRunningStats computes scalar mean and stddev",
    "[infra][stats][running_stats]")
{
    ScalarRunningStats stats;
    static_assert(
        std::same_as<decltype(stats.result()), ScalarRunningStatsResult>);

    SECTION("single sample")
    {
        stats.update(3.0);

        const auto result = stats.result();
        REQUIRE(result.mean == Approx(3.0));
        REQUIRE(result.stddev == Approx(0.0));
    }

    SECTION("multiple samples")
    {
        stats.update(2.0);
        stats.update(4.0);

        const auto result = stats.result();
        REQUIRE(result.mean == Approx(3.0));
        REQUIRE(result.stddev == Approx(std::sqrt(2.0)));
    }
}

TEST_CASE(
    "VectorRunningStats computes vector mean and stddev",
    "[infra][stats][running_stats]")
{
    SECTION("requires a positive size")
    {
        REQUIRE_THROWS_AS(VectorRunningStats{0}, GelexException);
        REQUIRE_THROWS_AS(VectorRunningStats{-1}, GelexException);
    }

    VectorRunningStats stats{2};
    static_assert(
        std::same_as<decltype(stats.result()), VectorRunningStatsResult>);

    SECTION("single sample")
    {
        stats.update(Eigen::VectorXd{{1.0, 3.0}});

        const auto result = stats.result();
        REQUIRE(result.mean.isApprox(Eigen::VectorXd{{1.0, 3.0}}));
        REQUIRE(result.stddev.isApprox(Eigen::VectorXd::Zero(2)));
    }

    SECTION("multiple samples")
    {
        stats.update(Eigen::VectorXd{{1.0, 3.0}});
        stats.update(Eigen::VectorXd{{3.0, 7.0}});

        const auto result = stats.result();
        REQUIRE(result.mean.isApprox(Eigen::VectorXd{{2.0, 5.0}}));
        REQUIRE(result.stddev.isApprox(
            Eigen::VectorXd{{std::sqrt(2.0), std::sqrt(8.0)}}));
    }
}

TEST_CASE(
    "CategoricalFrequency computes category probabilities",
    "[infra][stats][running_stats]")
{
    CategoricalFrequency frequency{3, 4};

    SECTION("no update gives zero probabilities")
    {
        const auto result = std::move(frequency).take_probabilities();
        REQUIRE(result.value.isApprox(Eigen::MatrixXd::Zero(3, 4)));
    }

    SECTION("counts normalized by update count")
    {
        frequency.update(Eigen::VectorXi{{0, 1, 3}});
        frequency.update(Eigen::VectorXi{{1, 1, 0}});
        frequency.update(Eigen::VectorXi{{1, 2, 3}});

        const Eigen::MatrixXd expected{
            {1.0 / 3.0, 2.0 / 3.0, 0.0, 0.0},
            {0.0, 2.0 / 3.0, 1.0 / 3.0, 0.0},
            {1.0 / 3.0, 0.0, 0.0, 2.0 / 3.0}};
        const auto result = std::move(frequency).take_probabilities();
        REQUIRE(result.value.isApprox(expected));
    }
}

}  // namespace gelex
