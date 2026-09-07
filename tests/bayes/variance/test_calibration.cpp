// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <limits>

#include "gelex/bayes/variance/calibration.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

using Catch::Approx;

TEST_CASE(
    "mean-calibrated variance parameter owns its initial value and prior",
    "[bayes][variance][calibration]")
{
    constexpr double target = 2.0;
    const auto parameter
        = gelex::detail::make_mean_calibrated_variance_parameter(target);
    const auto prior = parameter.prior.scaled_inv_chi2_parameters();

    REQUIRE(parameter.initial == target);
    REQUIRE(prior.degrees_of_freedom() == 4.0);
    REQUIRE(prior.scale() == Approx(0.5 * target));
}

TEST_CASE(
    "mean-calibrated variance parameter rejects invalid targets",
    "[bayes][variance][calibration]")
{
    SECTION("non-positive")
    {
        REQUIRE_THROWS_AS(
            gelex::detail::make_mean_calibrated_variance_parameter(0.0),
            gelex::GelexException);
    }

    SECTION("non-finite")
    {
        REQUIRE_THROWS_AS(
            gelex::detail::make_mean_calibrated_variance_parameter(
                std::numeric_limits<double>::infinity()),
            gelex::GelexException);
    }
}

TEST_CASE(
    "marker variance calibration uses owned per-mode values",
    "[bayes][variance][calibration]")
{
    std::array<double, 2> values{2.0, 3.0};
    const gelex::MarkerVarianceCalibrator calibrator{values};
    values.fill(100.0);
    REQUIRE(calibrator.calibrate(gelex::GeneticMode::A, 0.5).initial == 4.0);
    REQUIRE(calibrator.calibrate(gelex::GeneticMode::D, 0.25).initial == 12.0);
}

TEST_CASE(
    "marker variance calibration rejects invalid numeric inputs",
    "[bayes][variance][calibration]")
{
    const std::array invalid_values{
        -1.0,
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::quiet_NaN()};
    const gelex::MarkerVarianceCalibrator calibrator{{2.0, 0.0}};
    for (const double invalid : invalid_values)
    {
        REQUIRE_THROWS_AS(
            (gelex::MarkerVarianceCalibrator{{invalid, 1.0}}),
            gelex::GelexException);
        REQUIRE_THROWS_AS(
            calibrator.calibrate(gelex::GeneticMode::A, invalid),
            gelex::GelexException);
    }
    REQUIRE_THROWS_AS(
        calibrator.calibrate(gelex::GeneticMode::A, 0.0),
        gelex::GelexException);
    REQUIRE_THROWS_AS(
        calibrator.calibrate(gelex::GeneticMode::D, 1.0),
        gelex::GelexException);
}
