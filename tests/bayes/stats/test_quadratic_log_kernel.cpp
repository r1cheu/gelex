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
#include <numbers>

#include "gelex/bayes/stats/quadratic_log_kernel.h"

using Catch::Approx;
using gelex::make_coefficient_likelihood;
using gelex::make_normal_prior;
using gelex::QuadraticLogKernel;

TEST_CASE(
    "QuadraticLogKernel exposes its normalized normal parameters",
    "[bayes][stats][quadratic_log_kernel]")
{
    const QuadraticLogKernel kernel{3.0, 6.0, -4.0};
    const auto parameters = kernel.normal_parameters();

    const Eigen::Vector2d actual{parameters.mean(), parameters.stddev()};
    const Eigen::Vector2d expected{2.0, std::sqrt(1.0 / 3.0)};
    REQUIRE(actual.isApprox(expected));
}

TEST_CASE(
    "QuadraticLogKernel combines likelihood and normalized prior",
    "[bayes][stats][quadratic_log_kernel]")
{
    const QuadraticLogKernel likelihood{2.5, -1.5, 0.0};
    const QuadraticLogKernel prior{
        0.25, 0.0, -0.5 * std::log(8.0 * std::numbers::pi)};
    const auto posterior = likelihood + prior;
    const auto parameters = posterior.normal_parameters();

    const Eigen::Vector2d actual{parameters.mean(), parameters.stddev()};
    const Eigen::Vector2d expected{-6.0 / 11.0, 2.0 / std::sqrt(11.0)};
    REQUIRE(actual.isApprox(expected));

    const double expected_log_integral = -0.5 * (std::log(11.0) - (9.0 / 11.0));
    REQUIRE(posterior.log_integral() == Approx(expected_log_integral));
}

TEST_CASE(
    "make_normal_prior produces a normalized zero-mean kernel",
    "[bayes][stats][quadratic_log_kernel]")
{
    const auto prior = make_normal_prior(4.0);
    const auto parameters = prior.normal_parameters();

    const Eigen::Vector2d actual{parameters.mean(), parameters.stddev()};
    const Eigen::Vector2d expected{0.0, 2.0};
    REQUIRE(actual.isApprox(expected));
    REQUIRE(prior.log_integral() == Approx(0.0).margin(1e-12));
}

TEST_CASE(
    "make_coefficient_likelihood applies residual precision",
    "[bayes][stats][quadratic_log_kernel]")
{
    const auto likelihood = make_coefficient_likelihood(8.0, 6.0, 2.0);
    const auto parameters = likelihood.normal_parameters();

    const Eigen::Vector2d actual{parameters.mean(), parameters.stddev()};
    const Eigen::Vector2d expected{0.75, 0.5};
    REQUIRE(actual.isApprox(expected));
}

TEST_CASE(
    "QuadraticLogKernel log integral avoids avoidable intermediate overflow",
    "[bayes][stats][quadratic_log_kernel]")
{
    const QuadraticLogKernel large_linear{1e200, 1e200, 0.0};
    const QuadraticLogKernel small_quadratic{1e-310, 0.0, 0.0};

    REQUIRE(std::isfinite(large_linear.log_integral()));
    REQUIRE(std::isfinite(small_quadratic.log_integral()));
}
