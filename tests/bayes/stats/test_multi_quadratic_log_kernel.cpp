// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>

#include "gelex/bayes/stats/multi_quadratic_log_kernel.h"

using Catch::Approx;
using gelex::make_multi_normal_prior;
using gelex::MultiQuadraticLogKernel;

TEST_CASE(
    "MultiQuadraticLogKernel combines a likelihood and normal prior",
    "[bayes][stats][multi_quadratic_log_kernel]")
{
    const MultiQuadraticLogKernel likelihood{
        Eigen::Matrix2d{{1.0, 2.0}, {2.0, 4.0}},
        Eigen::Vector2d{{3.0, 6.0}},
        0.0};
    const auto prior
        = make_multi_normal_prior(Eigen::Matrix2d{{2.0, 1.0}, {1.0, 2.0}});

    const auto posterior = likelihood + prior;

    const Eigen::Matrix2d expected_quadratic{
        {5.0 / 3.0, 5.0 / 3.0}, {5.0 / 3.0, 14.0 / 3.0}};
    REQUIRE(posterior.quadratic().isApprox(expected_quadratic));
    REQUIRE(posterior.linear().isApprox(Eigen::Vector2d{{3.0, 6.0}}));
    REQUIRE(
        posterior.constant()
        == Approx(-std::log(2.0 * std::numbers::pi) - (0.5 * std::log(3.0))));
}
