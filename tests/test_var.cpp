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
#include <catch2/catch_test_macros.hpp>

#include "gelex/infra/stats/detail/var.h"

namespace gelex
{

TEST_CASE("vecvar computes vector variance", "[stats][var]")
{
    const Eigen::VectorXd values{{1.0, 2.0, 4.0, 7.0}};

    REQUIRE(stats::detail::vecvar(values) == 7.0);
    REQUIRE(stats::detail::vecvar(values, 0) == 5.25);
}

TEST_CASE("matvar computes axis-wise matrix variance", "[stats][var]")
{
    const Eigen::MatrixXd values{
        {1.0, 2.0},
        {3.0, 4.0},
        {5.0, 8.0},
    };
    const Eigen::RowVectorXd expected_colwise{{4.0, 28.0 / 3.0}};
    const Eigen::VectorXd expected_rowwise{{0.5, 0.5, 4.5}};

    REQUIRE(stats::detail::matvar<0>(values).isApprox(expected_colwise));
    REQUIRE(stats::detail::matvar<1>(values).isApprox(expected_rowwise));
}

}  // namespace gelex
