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
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <span>

#include "gelex/bayes/draw.h"
#include "gelex/infra/stats/result.h"
#include "gelex/io/binary_reader.h"
#include "gelex/io/binary_writer.h"

#include "file_fixture.h"

namespace gelex
{

using Catch::Approx;
TEST_CASE(
    "ScalarDraw writes retained samples and computes statistics",
    "[bayes][draw]")
{
    test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "scalar.draws";
    ScalarRunningStatsResult result;

    {
        BinaryWriter writer(path.string());
        ScalarDraw draw{writer.reserve<double>("scalar", BinaryShape{1, 2})};

        draw.append(2.0);
        draw.append(4.0);
        writer.close();
        result = draw.result();
    }

    REQUIRE(result.mean == Approx(3.0));
    REQUIRE(result.stddev == Approx(std::sqrt(2.0)));

    const BinaryReader reader(path.string());
    REQUIRE(
        reader.to_map<double>("scalar").isApprox(Eigen::MatrixXd{{2.0, 4.0}}));
}

TEST_CASE(
    "VectorDraw writes retained samples and computes statistics",
    "[bayes][draw]")
{
    test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "vector.draws";
    VectorRunningStatsResult result;

    {
        BinaryWriter writer(path.string());
        VectorDraw draw{writer.reserve<float>("vector", BinaryShape{2, 2})};

        const std::array first{1.0, 3.0};
        draw.append(std::span{first});
        draw.append(Eigen::VectorXd{{3.0, 7.0}});
        writer.close();
        result = draw.result();
    }

    REQUIRE(result.mean.isApprox(Eigen::VectorXd{{2.0, 5.0}}));
    REQUIRE(result.stddev.isApprox(
        Eigen::VectorXd{{std::sqrt(2.0), std::sqrt(8.0)}}));

    const BinaryReader reader(path.string());
    REQUIRE(reader.to_map<float>("vector").isApprox(
        Eigen::MatrixXf{{1.0, 3.0}, {3.0, 7.0}}));
}

TEST_CASE(
    "CategoryDraw writes retained samples and computes probabilities",
    "[bayes][draw]")
{
    test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "category.draws";
    CategoryRunningStatsResult result;

    {
        BinaryWriter writer(path.string());
        CategoryDraw<3> draw{
            writer.reserve<std::uint8_t>("category", BinaryShape{2, 2})};

        const std::array<std::uint8_t, 2> first{0, 1};
        draw.append(std::span{first});
        draw.append(Eigen::VectorX<std::uint8_t>{{2, 1}});
        writer.close();
        result = draw.result();
    }

    REQUIRE(result.probabilities.isApprox(
        Eigen::MatrixXd{{0.5, 0.0, 0.5}, {0.0, 1.0, 0.0}}));

    const BinaryReader reader(path.string());
    REQUIRE(reader.to_map<std::uint8_t>("category")
                .isApprox(
                    Eigen::Matrix<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic>{
                        {0, 2}, {1, 1}}));
}

}  // namespace gelex
