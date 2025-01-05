#include <armadillo>
#include <catch2/catch_test_macros.hpp>
#include "chenx/data/encode.h"
#include "chenx/data/impute.h"

TEST_CASE("impute method work correct", "[vector]")
{
    // For each section, vector v is anew:

    arma::dmat X = {{NAN, 2.0, 3.0}, {4.0, NAN, 6.0}, {7.0, 8.0, NAN}};

    SECTION("mean imputation")
    {
        arma::dvec mean = chenx::MeanImpute(X);
        arma::dmat expected
            = {{5.5, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 4.5}};
        REQUIRE(arma::approx_equal(X, expected, "absdiff", 1e-10));
    }
    SECTION("median even imputation")
    {
        arma::dvec median = chenx::MedianImpute(X);
        arma::dmat expected
            = {{5.5, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 4.5}};
        REQUIRE(arma::approx_equal(X, expected, "absdiff", 1e-10));
    }
    SECTION("median odd imputation")
    {
        arma::dmat X_odd = {
            {NAN, 2.0, 3.0}, {4.0, NAN, 6.0}, {7.0, 8.0, NAN}, {1.0, 2.0, 3.0}};
        arma::dvec median = chenx::MedianImpute(X_odd);
        arma::dmat expected = {
            {4.0, 2.0, 3.0}, {4.0, 2.0, 6.0}, {7.0, 8.0, 3.0}, {1.0, 2.0, 3.0}};

        REQUIRE(arma::approx_equal(X_odd, expected, "absdiff", 1e-10));
    }
    SECTION("value imputation")
    {
        chenx::ValueImpute(X, {3.0, 3.0, 3.0});
        arma::dmat expected
            = {{3.0, 2.0, 3.0}, {4.0, 3.0, 6.0}, {7.0, 8.0, 3.0}};
        REQUIRE(arma::approx_equal(X, expected, "absdiff", 1e-10));
    }
}

TEST_CASE("EncodeTest Hybird", "[encode]")
{
    arma::dmat X = {
        {1, 0, 2, 2}, {1, 2, 2, 1}, {2, 2, 2, 2}, {2, 2, 2, 1}, {1, 0, 2, 2}};
    arma::dmat guide = {{0, 0, 0, 2}, {1, 1.5, 2, 2.5}};
    chenx::HybridEncode(X, guide);

    arma::dmat Y
        = {{1, 0, 2, 0},
           {1, 2, 2, 2.5},
           {2, 2, 2, 0},
           {2, 2, 2, 2.5},
           {1, 0, 2, 0}};
    REQUIRE(arma::approx_equal(X, Y, "absdiff", 1e-10));
}

TEST_CASE("ComputeHybirdValue Tests", "[hybird]")
{
    arma::dvec phenotype = {1.0, 2.0, 3.0, 4.0};

    SECTION("Basic")
    {
        arma::dmat X = {
            {0, 1, 2},
            {1, 0, 2},
            {2, 1, 0},
            {1, 2, 1},
        };
        arma::dmat result = chenx::ComputeHybirdValue(X, phenotype);
        arma::dmat expected = {{0, 0, 2}, {2, 0, 10 / 3.0}};
        REQUIRE(arma::approx_equal(result, expected, "absdiff", 1e-10));
    }

    SECTION("Not All [0 1 2] exist")
    {
        arma::dmat X = {
            {0, 1, 2},
            {1, 0, 2},
            {2, 1, 0},
            {1, 2, 0},
        };
        arma::dmat result = chenx::ComputeHybirdValue(X, phenotype);
        arma::dmat expected = {{0, 0, 0}, {2, 0, 1}};
        REQUIRE(arma::approx_equal(result, expected, "absdiff", 1e-10));
    }

    SECTION("NaNHandling")
    {
        arma::dmat X = {
            {0, 1, 2},
            {1, 0, arma::datum::nan},
            {2, 1, 0},
            {1, 2, 1},
        };
        arma::dmat result = chenx::ComputeHybirdValue(X, phenotype);
        arma::dmat expected = {{0, 0, 2}, {2, 0, 3}};
        REQUIRE(arma::approx_equal(result, expected, "absdiff", 1e-10));
    }
}
