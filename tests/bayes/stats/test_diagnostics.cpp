// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <random>

#include "gelex/bayes/stats/diagnostics.h"
#include "gelex/exception.h"

extern "C" const char* __lsan_default_suppressions()
{
    return "leak:___kmp_allocate_align\n"
           "leak:libomp\n"
           "leak:libiomp5";
}

using Catch::Matchers::WithinAbs;

TEST_CASE("diagnose chain", "[bayes][stats][diagnostics]")
{
    std::mt19937_64 rng(7);
    std::normal_distribution<double> dist(1.0, 2.0);
    const Eigen::VectorXd draws
        = Eigen::VectorXd::NullaryExpr(200, [&]() { return dist(rng); });

    const auto result = gelex::diagnose_chain(draws, 0.9);

    SECTION("descriptive statistics")
    {
        Eigen::VectorXd sorted = draws;
        std::sort(sorted.begin(), sorted.end());
        const double sd
            = std::sqrt((draws.array() - draws.mean()).square().sum() / 199.0);
        REQUIRE_THAT(result.mean, WithinAbs(draws.mean(), 1e-12));
        REQUIRE_THAT(result.sd, WithinAbs(sd, 1e-12));
        REQUIRE_THAT(
            result.median, WithinAbs((sorted(99) + sorted(100)) / 2.0, 1e-12));
        REQUIRE_THAT(result.mcse, WithinAbs(sd / std::sqrt(result.ess), 1e-12));
        REQUIRE(result.hpdi_lower < result.median);
        REQUIRE(result.median < result.hpdi_upper);
    }

    SECTION("hpdi of a full probability mass is the sample range")
    {
        const auto full = gelex::diagnose_chain(draws, 1.0);
        REQUIRE(full.hpdi_lower == draws.minCoeff());
        REQUIRE(full.hpdi_upper == draws.maxCoeff());
    }

    SECTION("hpdi is the narrowest interval")
    {
        std::exponential_distribution<double> exponential;
        const Eigen::VectorXd skewed = Eigen::VectorXd::NullaryExpr(
            20000, [&]() { return exponential(rng); });
        const auto tail = gelex::diagnose_chain(skewed, 0.2);
        REQUIRE_THAT(tail.hpdi_lower, WithinAbs(0.0, 0.01));
        REQUIRE_THAT(tail.hpdi_upper, WithinAbs(0.22, 0.01));
    }

    SECTION("effective sample size of a linear trend over chains")
    {
        const Eigen::VectorXd trend = Eigen::VectorXd::LinSpaced(1000, 0, 999);
        const Eigen::MatrixXd chains = trend.reshaped(10, 100);
        REQUIRE_THAT(gelex::diagnose_chain(chains).ess, WithinAbs(52.64, 0.2));
    }

    SECTION("independent draws have ess near the draw count")
    {
        REQUIRE(result.ess > 150.0);
        REQUIRE(result.ess < 260.0);
    }

    SECTION("split rhat of a stationary chain is near one")
    {
        REQUIRE_THAT(result.split_rhat, WithinAbs(1.0, 0.05));
    }

    SECTION("split rhat detects a shifted chain")
    {
        Eigen::MatrixXd chains(200, 2);
        chains.col(0) = draws;
        chains.col(1) = draws.array() + 10.0;
        REQUIRE(gelex::diagnose_chain(chains).split_rhat > 2.0);
    }

    SECTION("accepts any vector orientation and column layout")
    {
        Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(3, 200);
        mat.row(1) = draws.transpose();
        const Eigen::RowVectorXd row = draws.transpose();

        const auto from_row = gelex::diagnose_chain(row, 0.9);
        const auto from_mat_row = gelex::diagnose_chain(mat.row(1), 0.9);
        const auto from_mat_col
            = gelex::diagnose_chain(mat.transpose().col(1), 0.9);
        for (const auto& other : {from_row, from_mat_row, from_mat_col})
        {
            REQUIRE(other.mean == result.mean);
            REQUIRE(other.ess == result.ess);
            REQUIRE(other.split_rhat == result.split_rhat);
            REQUIRE(other.hpdi_upper == result.hpdi_upper);
        }
    }

    SECTION("constant chain yields NaN diagnostics")
    {
        const Eigen::VectorXd constant = Eigen::VectorXd::Constant(10, 3.0);
        const auto flat = gelex::diagnose_chain(constant);
        REQUIRE(flat.mean == 3.0);
        REQUIRE(flat.sd == 0.0);
        REQUIRE(flat.hpdi_lower == 3.0);
        REQUIRE(flat.hpdi_upper == 3.0);
        REQUIRE(std::isnan(flat.ess));
        REQUIRE(std::isnan(flat.mcse));
        REQUIRE(std::isnan(flat.split_rhat));
    }

    SECTION("rejects invalid input")
    {
        REQUIRE_THROWS_AS(
            gelex::diagnose_chain(Eigen::VectorXd{{1.0, 2.0, 3.0}}),
            gelex::GelexException);
        REQUIRE_THROWS_AS(
            gelex::diagnose_chain(Eigen::MatrixXd(10, 0)),
            gelex::GelexException);
        REQUIRE_THROWS_AS(
            gelex::diagnose_chain(draws, 0.0), gelex::GelexException);
        REQUIRE_THROWS_AS(
            gelex::diagnose_chain(draws, 1.5), gelex::GelexException);
    }
}
