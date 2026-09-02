// Copyright 2026 RuLei Chen
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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
