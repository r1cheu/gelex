// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <span>

#include "gelex/bayes/stats/detail/running_stats.h"
#include "gelex/bayes/stats/result.h"
#include "gelex/exception.h"

namespace gelex
{

using Catch::Approx;
using gelex::detail::CategoryRunningStats;
using gelex::detail::ScalarRunningStats;
using gelex::detail::VectorRunningStats;

TEST_CASE(
    "ScalarRunningStats computes scalar mean and stddev",
    "[bayes][stats][running_stats]")
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
    "[bayes][stats][running_stats]")
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
        REQUIRE(stats.mean_square().isApprox(Eigen::VectorXd{{5.0, 29.0}}));
    }
}

TEST_CASE(
    "CategoryRunningStats computes per-item category probabilities",
    "[bayes][stats][running_stats]")
{
    SECTION("requires a positive size")
    {
        REQUIRE_THROWS_AS(CategoryRunningStats<3>{0}, GelexException);
        REQUIRE_THROWS_AS(CategoryRunningStats<3>{-1}, GelexException);
    }

    CategoryRunningStats<3> stats{2};
    static_assert(
        std::same_as<decltype(stats.result()), CategoryRunningStatsResult>);

    const std::array<std::uint8_t, 2> first{0, 1};
    stats.update(first);

    SECTION("single sample")
    {
        const auto result = stats.result();
        const Eigen::MatrixXd expected{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
        REQUIRE(result.probabilities.isApprox(expected));
    }

    SECTION("multiple samples")
    {
        const std::array<std::uint8_t, 2> second{2, 1};
        stats.update(std::span{second});

        const auto result = stats.result();
        const Eigen::MatrixXd expected{{0.5, 0.0, 0.5}, {0.0, 1.0, 0.0}};
        REQUIRE(result.probabilities.isApprox(expected));
    }
}

}  // namespace gelex
