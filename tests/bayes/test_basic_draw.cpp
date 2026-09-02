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
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <cstdint>
#include <span>
#include <utility>

#include "gelex/bayes/basic_draw.h"
#include "gelex/bayes/stats/result.h"
#include "gelex/io/binary_reader.h"
#include "gelex/io/binary_writer.h"

#include "file_fixture.h"

namespace gelex
{

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;

TEST_CASE("Basic draws expose their payload identifiers", "[bayes][draw]")
{
    test::FileFixture fixture;
    BinaryWriter writer(
        (fixture.get_test_dir() / "draw_identifiers.draws").string());

    ScalarDraw scalar{writer.reserve<double>("scalar/path", BinaryShape{1, 0})};
    VectorDraw vector{writer.reserve<float>("vector/path", BinaryShape{2, 0})};
    CategoryDraw<3> category{
        writer.reserve<std::uint8_t>("category/path", BinaryShape{2, 0})};

    ScalarDraw moved_scalar{std::move(scalar)};
    VectorDraw moved_vector{std::move(vector)};
    CategoryDraw<3> moved_category{std::move(category)};

    REQUIRE(moved_scalar.identifier() == "scalar/path");
    REQUIRE(moved_vector.identifier() == "vector/path");
    REQUIRE(moved_category.identifier() == "category/path");
}

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

TEST_CASE(
    "Draws reject statistics requested before any sample",
    "[bayes][draw]")
{
    test::FileFixture fixture;
    BinaryWriter writer(
        (fixture.get_test_dir() / "empty_draws.draws").string());

    const ScalarDraw scalar{
        writer.reserve<double>("scalar/path", BinaryShape{1, 1})};
    const VectorDraw vector{
        writer.reserve<float>("vector/path", BinaryShape{2, 1})};
    const CategoryDraw<3> category{
        writer.reserve<std::uint8_t>("category/path", BinaryShape{2, 1})};

    REQUIRE_THROWS_WITH(
        scalar.result(),
        ContainsSubstring("posterior 'scalar/path' has no recorded draws"));
    REQUIRE_THROWS_WITH(
        vector.result(),
        ContainsSubstring("posterior 'vector/path' has no recorded draws"));
    REQUIRE_THROWS_WITH(
        category.result(),
        ContainsSubstring("posterior 'category/path' has no recorded draws"));
}

}  // namespace gelex
