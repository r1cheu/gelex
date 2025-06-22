#include "gelex/estimator/bayes/diagnostics.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <random>
#include "armadillo"

using arma::approx_equal;
using arma::dcube;
using arma::dmat;
using arma::zeros;

TEST_CASE("gelman rubin", "[diagnostics]")
{
    dcube x = zeros(1, 10, 2);
    x.slice(0) = arma::linspace(0, 9, 10).t();
    x.slice(1) = arma::linspace(0, 9, 10).t() + 1;

    double rhat = as_scalar(gelex::gelman_rubin(x));
    REQUIRE(rhat == Catch::Approx(0.98).margin(0.01));
}

TEST_CASE("split gelman rubin", "[diagnostics]")
{
    std::mt19937_64 rng(42);
    std::normal_distribution<double> dist(0, 1);
    dcube x(1, 10, 2);
    x.imbue([&]() { return dist(rng); });
    dcube a = x;
    x.reshape(5, 2, 2).reshape(1, 5, 4);

    auto r_hat1 = gelex::gelman_rubin(x);
    auto r_hat2 = gelex::split_gelman_rubin(a);

    REQUIRE(approx_equal(r_hat1, r_hat2, "absdiff", 1e-7));
}

TEST_CASE("autocorrelation", "[diagnostics]")
{
    dcube x(1, 10, 1);
    x.slice(0) = arma::linspace(0, 9, 10).t();

    dcube random_x(1, 20000, 1);
    std::mt19937_64 rng(42);
    std::normal_distribution<double> dist(0, 1);
    random_x.imbue([&]() { return dist(rng); });

    SECTION("autocorrelation with unbiased estimator")
    {
        dcube expected(1, 10, 1);
        expected.slice(0).row(0)
            = {1, 0.78, 0.52, 0.21, -0.13, -0.52, -0.94, -1.4, -1.91, -2.45};
        REQUIRE(approx_equal(
            gelex::autocorrelation(x, false), expected, "absdiff", 0.01));
    }

    SECTION("autocorrelation with biased estimator")
    {
        dcube expected(1, 10, 1);
        expected.slice(0).row(0)
            = {1, 0.78, 0.52, 0.21, -0.13, -0.52, -0.94, -1.4, -1.91, -2.45};

        expected.slice(0).row(0) %= arma::regspace(10, -1, 1).t();
        expected /= 10;
        REQUIRE(approx_equal(
            gelex::autocorrelation(x, true), expected, "absdiff", 0.01));
    }

    SECTION("autocorrelation with random data and unbiased estimator")
    {
        dcube result = gelex::autocorrelation(random_x, false);
        arma::rowvec ac
            = result.slice(0).cols(result.n_cols - 100, result.n_cols - 1);
        REQUIRE(arma::any(arma::abs(ac) > 0.1));
    }

    SECTION("autocorrelation with random data and biased estimator")
    {
        dcube result = gelex::autocorrelation(random_x, true);
        arma::rowvec ac
            = result.slice(0).cols(result.n_cols - 100, result.n_cols - 1);
        REQUIRE(approx_equal(
            arma::abs(ac), arma::zeros(arma::size(ac)), "absdiff", 0.01));
    }
}

TEST_CASE("autovariance", "[diagnostics]")
{
    dcube x(1, 10, 1);
    x.slice(0) = arma::linspace(0, 9, 10).t();
    auto actual = gelex::autocovariance(x, false);
    x.slice(0)
        = {8.25, 6.42, 4.25, 1.75, -1.08, -4.25, -7.75, -11.58, -15.75, -20.25};

    REQUIRE(approx_equal(actual, x, "absdiff", 0.01));
}

TEST_CASE("effective sample size", "[diagnostics]")
{
    dcube x(1, 1000, 1);
    x.slice(0) = arma::linspace(0, 999, 1000).t();
    x.reshape(1, 10, 100);

    auto result = gelex::effect_sample_size(x, false);
    REQUIRE(approx_equal(result, arma::dvec{52.64}, "absdiff", 0.01));
}

TEST_CASE("hpdi", "[diagnostics]")
{
    dcube x(2, 20000, 2);
    std::mt19937_64 rng(42);
    std::exponential_distribution<double> dist;
    x.imbue([&]() { return dist(rng); });

    dmat expected_exp = {{0.0, 0.22}, {0.0, 0.22}};
    dmat result = gelex::hpdi(x, 0.2);
    REQUIRE(approx_equal(result, expected_exp, "absdiff", 0.01));
}
