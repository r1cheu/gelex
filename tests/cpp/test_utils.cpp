#include <armadillo>
#include <catch2/catch_test_macros.hpp>
#include "gelex/utils/utils.h"

using namespace arma;
using namespace gelex;

TEST_CASE("centralize works correctly", "[utils]")
{
    dmat x = {{1, 4}, {2, 5}, {3, 6}};
    dvec expected_means = {2, 5};
    dmat expected_x = {{-1, -1}, {0, 0}, {1, 1}};

    dvec actual_means = centralize(x);

    REQUIRE(approx_equal(actual_means, expected_means, "absdiff", 1e-7));
    REQUIRE(approx_equal(x, expected_x, "absdiff", 1e-7));
}

TEST_CASE("standradize works correctly", "[utils]")
{
    dmat x = {{1, 4}, {2, 5}, {3, 6}};
    dvec expected_means = {2, 5};
    dvec expected_stddevs = {1.0, 1.0};  // After standardization
    dmat expected_x = {{-1, -1}, {0, 0}, {1, 1}};

    auto [actual_means, actual_stddevs] = standradize(x);

    REQUIRE(approx_equal(actual_means, expected_means, "absdiff", 1e-7));
    REQUIRE(approx_equal(actual_stddevs, expected_stddevs, "absdiff", 1e-7));
    REQUIRE(approx_equal(x, expected_x, "absdiff", 1e-7));
}

TEST_CASE("standradize handles constant columns", "[utils]")
{
    dmat x = {{1, 5}, {1, 5}, {1, 5}};
    dvec expected_means = {1, 5};
    dvec expected_stddevs = {0, 0};
    dmat expected_x(3, 2, fill::zeros);  // All zeros

    auto [means, stddevs] = standradize(x);

    REQUIRE(approx_equal(means, expected_means, "absdiff", 1e-7));
    REQUIRE(approx_equal(stddevs, expected_stddevs, "absdiff", 1e-7));
    REQUIRE(approx_equal(x, expected_x, "absdiff", 1e-7));
}
