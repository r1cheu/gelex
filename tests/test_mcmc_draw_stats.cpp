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

#include <cmath>
#include <utility>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/algo/mcmc/detail/records.h"

namespace gelex
{

TEST_CASE("ScalarDrawStats accumulates scalar draws", "[mcmc][draw_stats]")
{
    mcmc::detail::ScalarDrawStats draws;

    draws.store(2.0);
    draws.store(4.0);

    const auto result = draws.stats();
    REQUIRE(result.mean.isApprox(Eigen::VectorXd{{3.0}}));
    REQUIRE(result.stddev.isApprox(Eigen::VectorXd{{std::sqrt(2.0)}}));
}

TEST_CASE("VectorDrawStats accumulates vector draws", "[mcmc][draw_stats]")
{
    mcmc::detail::VectorDrawStats draws;

    const Eigen::VectorXd first{{1.0, 3.0}};
    const Eigen::VectorXd second{{3.0, 7.0}};
    draws.store(first);
    draws.store(second);

    const auto result = draws.stats();
    REQUIRE(result.mean.isApprox(Eigen::VectorXd{{2.0, 5.0}}));
    REQUIRE(
        result.stddev.isApprox(
            Eigen::VectorXd{{std::sqrt(2.0), std::sqrt(8.0)}}));
}

TEST_CASE(
    "CategoricalDrawStats accumulates categorical draws",
    "[mcmc][draw_stats]")
{
    mcmc::detail::CategoricalDrawStats draws{3, 4};

    draws.store(Eigen::VectorXi{{1, 0, 3}});
    draws.store(Eigen::VectorXi{{3, 0, 0}});

    const Eigen::MatrixXd expected_probabilities{
        {0.0, 0.5, 0.0, 0.5},
        {1.0, 0.0, 0.0, 0.0},
        {0.5, 0.0, 0.0, 0.5}};
    const auto result = std::move(draws).take_probabilities();
    REQUIRE(result.value.isApprox(expected_probabilities));
}

}  // namespace gelex
