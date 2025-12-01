
#include <iostream>
#include <random>

#include <Eigen/Dense>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../src/estimator/bayes/diagnostics.h"
#include "gelex/types/mcmc_samples.h"

using Catch::Matchers::WithinAbs;
using Eigen::RowVectorXd;

TEST_CASE("gelman rubin", "[diagnostics]")
{
    gelex::Samples samples(2);
    samples[0] = Eigen::RowVectorXd::LinSpaced(10, 0, 9);
    samples[1] = Eigen::RowVectorXd::LinSpaced(10, 0, 9).array() + 1;

    double rhat = gelex::gelman_rubin(samples)(0);
    REQUIRE_THAT(rhat, WithinAbs(0.98, 0.01));
}

TEST_CASE("split gelman rubin", "[diagnostics]")
{
    std::mt19937_64 rng(42);
    std::normal_distribution<double> dist(0, 1);

    gelex::Samples samples(2);
    samples[0]
        = Eigen::RowVectorXd::NullaryExpr(10, [&]() { return dist(rng); });
    samples[1]
        = Eigen::RowVectorXd::NullaryExpr(10, [&]() { return dist(rng); });
    gelex::Samples repeated_samples(4);

    for (Eigen::Index i = 0; i < 2; ++i)
    {
        repeated_samples[2 * i] = samples[i].leftCols(5);
        repeated_samples[(2 * i) + 1] = samples[i].rightCols(5);
    }

    auto r_hat1 = gelex::gelman_rubin(repeated_samples);
    auto r_hat2 = gelex::split_gelman_rubin(samples);

    REQUIRE(r_hat1.isApprox(r_hat2, 1e-7));
}

TEST_CASE("autocorrelation", "[diagnostics]")
{
    gelex::Samples x(1);
    x[0] = Eigen::RowVectorXd::LinSpaced(10, 0, 9);

    gelex::Samples random_x(1);
    std::mt19937_64 rng(42);
    std::normal_distribution<double> dist(0, 1);
    random_x[0]
        = Eigen::RowVectorXd::NullaryExpr(20000, [&]() { return dist(rng); });

    SECTION("autocorrelation with unbiased estimator")
    {
        Eigen::RowVectorXd expected{
            {1, 0.78, 0.52, 0.21, -0.13, -0.52, -0.94, -1.4, -1.91, -2.45}};

        auto result = gelex::autocorrelation(x, false);
        REQUIRE(result[0].isApprox(expected, 0.01));
    }

    SECTION("autocorrelation with biased estimator")
    {
        Eigen::RowVectorXd expected{
            {1, 0.78, 0.52, 0.21, -0.13, -0.52, -0.94, -1.4, -1.91, -2.45}};

        RowVectorXd weights = Eigen::VectorXd::LinSpaced(10, 10, 1);
        expected = expected.array() * weights.array() / 10.0;

        auto result = gelex::autocorrelation(x, true);
        REQUIRE(result[0].isApprox(expected, 0.01));
    }

    SECTION("autocorrelation with random data and unbiased estimator")
    {
        auto result = gelex::autocorrelation(random_x, false);
        Eigen::RowVectorXd ac = result[0].rightCols(100);
        REQUIRE(ac.array().abs().maxCoeff() > 0.1);
    }

    SECTION("autocorrelation with random data and biased estimator")
    {
        auto result = gelex::autocorrelation(random_x, true);
        Eigen::RowVectorXd ac = result[0].rightCols(100);

        REQUIRE(ac.array().abs().maxCoeff() < 0.01);
    }
}

TEST_CASE("autovariance", "[diagnostics]")
{
    gelex::Samples x;
    x.emplace_back(Eigen::RowVectorXd::LinSpaced(10, 0, 9));
    auto actual = gelex::autocovariance(x, false);

    Eigen::RowVectorXd expected{
        {8.25, 6.42, 4.25, 1.75, -1.08, -4.25, -7.75, -11.58, -15.75, -20.25}};

    std::cout << expected << "\n";
    std::cout << actual[0].row(0) << "\n";
    REQUIRE(actual[0].row(0).isApprox(expected, 0.01));
}

TEST_CASE("effective sample size", "[diagnostics]")
{
    gelex::Samples x;

    Eigen::RowVectorXd raw_x = Eigen::RowVectorXd::LinSpaced(1000, 0, 999);
    for (Eigen::Index i = 0; i < 100; ++i)
    {
        x.emplace_back(raw_x.segment(i * 10, 10));
    }

    auto result = gelex::effect_sample_size(x, false);
    REQUIRE(std::abs(result(0) - 52.64) < 0.01);
}

TEST_CASE("hpdi", "[diagnostics]")
{
    Eigen::VectorXd x(20000);
    std::mt19937_64 rng(42);
    std::exponential_distribution<double> dist;
    for (int i = 0; i < x.size(); ++i)
    {
        x(i) = dist(rng);
    }

    auto [left, right] = gelex::hpdi(x, 0.2);
    REQUIRE_THAT(left, WithinAbs(0.0, 0.01));
    REQUIRE_THAT(right, WithinAbs(0.22, 0.01));
}
