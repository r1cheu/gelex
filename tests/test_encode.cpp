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

#include <format>
#include <string>
#include <type_traits>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/data/dataframe/column.h"
#include "gelex/data/dataframe/constants.h"
#include "gelex/data/dataframe/encode.h"
#include "gelex/exception.h"

using gelex::df::check_levels;
using gelex::df::collect_levels;
using gelex::df::Column;
using gelex::df::dummy_encode;
using gelex::df::encode;
using gelex::df::EncodedResult;
using gelex::df::kSeparator;
using gelex::df::LevelMismatch;
using gelex::df::one_hot_encode;

namespace
{

auto make_col(std::string_view name, std::vector<std::string> values) -> Column
{
    return Column(name, std::move(values));
}

}  // namespace

// ================================================================
// one_hot_encode
// ================================================================

TEST_CASE(
    "one_hot_encode three levels produces 3-column matrix",
    "[encode][dataframe]")
{
    auto col = make_col("group", {"A", "B", "C"});
    auto result = one_hot_encode(col);

    REQUIRE(result.data.rows() == 3);
    REQUIRE(result.data.cols() == 3);
    REQUIRE(result.name == "group");
    REQUIRE(result.level_names.size() == 3);
}

TEST_CASE("one_hot_encode levels are sorted", "[encode][dataframe]")
{
    auto col = make_col("g", {"B", "A", "B", "C"});
    auto result = one_hot_encode(col);

    REQUIRE(result.level_names[0] == std::format("g{}A", kSeparator));
    REQUIRE(result.level_names[1] == std::format("g{}B", kSeparator));
    REQUIRE(result.level_names[2] == std::format("g{}C", kSeparator));
}

TEST_CASE("one_hot_encode matrix has correct values", "[encode][dataframe]")
{
    // values: A, B, A, C -> sorted levels: A=col0, B=col1, C=col2
    auto col = make_col("g", {"A", "B", "A", "C"});
    auto result = one_hot_encode(col);

    Eigen::MatrixXd expected(4, 3);
    // clang-format off
    expected << 1, 0, 0,
                0, 1, 0,
                1, 0, 0,
                0, 0, 1;
    // clang-format on
    REQUIRE(result.data.isApprox(expected));
}

TEST_CASE("one_hot_encode single level throws", "[encode][dataframe]")
{
    auto col = make_col("x", {"only", "only", "only"});
    REQUIRE_THROWS_AS(one_hot_encode(col), gelex::InvalidInputException);
}

TEST_CASE(
    "one_hot_encode float scalar produces float matrix",
    "[encode][dataframe]")
{
    auto col = make_col("f", {"X", "Y"});
    auto result = one_hot_encode<float>(col);

    using MatType = decltype(result.data);
    static_assert(std::is_same_v<MatType::Scalar, float>);

    Eigen::MatrixXf expected(2, 2);
    // clang-format off
    expected << 1, 0,
                0, 1;
    // clang-format on
    REQUIRE(result.data.isApprox(expected));
}

// ================================================================
// dummy_encode
// ================================================================

TEST_CASE(
    "dummy_encode three levels produces 2-column matrix",
    "[encode][dataframe]")
{
    auto col = make_col("group", {"A", "B", "C"});
    auto result = dummy_encode(col);

    REQUIRE(result.data.rows() == 3);
    REQUIRE(result.data.cols() == 2);
    REQUIRE(result.name == "group");
    REQUIRE(result.level_names.size() == 2);
}

TEST_CASE(
    "dummy_encode drops first sorted level as reference",
    "[encode][dataframe]")
{
    // sorted: A, B, C — reference = A, active = {B, C}
    auto col = make_col("g", {"C", "B", "A"});
    auto result = dummy_encode(col);

    REQUIRE(result.level_names[0] == std::format("g{}B", kSeparator));
    REQUIRE(result.level_names[1] == std::format("g{}C", kSeparator));
}

TEST_CASE("dummy_encode matrix has correct values", "[encode][dataframe]")
{
    // sorted: A(ref), B=col0, C=col1
    auto col = make_col("g", {"A", "B", "C", "A"});
    auto result = dummy_encode(col);

    Eigen::MatrixXd expected(4, 2);
    // clang-format off
    expected << 0, 0,
                1, 0,
                0, 1,
                0, 0;
    // clang-format on
    REQUIRE(result.data.isApprox(expected));
}

TEST_CASE(
    "dummy_encode two levels produces 1-column matrix",
    "[encode][dataframe]")
{
    auto col = make_col("sex", {"M", "F", "M"});
    auto result = dummy_encode(col);

    REQUIRE(result.data.rows() == 3);
    REQUIRE(result.data.cols() == 1);
    REQUIRE(result.level_names[0] == std::format("sex{}M", kSeparator));
}

TEST_CASE("dummy_encode single level throws", "[encode][dataframe]")
{
    auto col = make_col("mono", {"X", "X", "X"});
    REQUIRE_THROWS_AS(dummy_encode(col), gelex::InvalidInputException);
}

// ================================================================
// encode
// ================================================================

TEST_CASE(
    "encode basic usage with user-specified levels",
    "[encode][dataframe]")
{
    auto col = make_col("g", {"A", "B", "C"});
    std::vector<std::string> levels = {"A", "B", "C"};
    auto result = encode(col, levels);

    REQUIRE(result.data.rows() == 3);
    REQUIRE(result.data.cols() == 3);
    REQUIRE(result.name == "g");

    Eigen::MatrixXd expected(3, 3);
    // clang-format off
    expected << 1, 0, 0,
                0, 1, 0,
                0, 0, 1;
    // clang-format on
    REQUIRE(result.data.isApprox(expected));
}

TEST_CASE(
    "encode preserves user-specified level order without sorting",
    "[encode][dataframe]")
{
    auto col = make_col("g", {"A", "B", "C"});
    std::vector<std::string> levels = {"C", "A", "B"};
    auto result = encode(col, levels);

    // C=col0, A=col1, B=col2
    REQUIRE(result.level_names[0] == std::format("g{}C", kSeparator));
    REQUIRE(result.level_names[1] == std::format("g{}A", kSeparator));
    REQUIRE(result.level_names[2] == std::format("g{}B", kSeparator));

    Eigen::MatrixXd expected(3, 3);
    // clang-format off
    expected << 0, 1, 0,
                0, 0, 1,
                1, 0, 0;
    // clang-format on
    REQUIRE(result.data.isApprox(expected));
}

TEST_CASE(
    "encode unknown values in data produce all-zero rows",
    "[encode][dataframe]")
{
    // "D" not in levels
    auto col = make_col("g", {"A", "D", "B"});
    std::vector<std::string> levels = {"A", "B"};
    auto result = encode(col, levels);

    Eigen::MatrixXd expected(3, 2);
    // clang-format off
    expected << 1, 0,
                0, 0,
                0, 1;
    // clang-format on
    REQUIRE(result.data.isApprox(expected));
}

TEST_CASE(
    "encode extra levels not in data produce all-zero columns",
    "[encode][dataframe]")
{
    // "Z" not in data
    auto col = make_col("g", {"A", "B"});
    std::vector<std::string> levels = {"A", "Z", "B"};
    auto result = encode(col, levels);

    REQUIRE(result.data.cols() == 3);

    Eigen::MatrixXd expected(2, 3);
    // clang-format off
    expected << 1, 0, 0,
                0, 0, 1;
    // clang-format on
    REQUIRE(result.data.isApprox(expected));
}

TEST_CASE("encode duplicate levels throws", "[encode][dataframe]")
{
    auto col = make_col("g", {"A", "B"});
    std::vector<std::string> levels = {"A", "B", "A"};
    REQUIRE_THROWS_AS(encode(col, levels), gelex::InvalidInputException);
}

TEST_CASE("encode empty levels throws", "[encode][dataframe]")
{
    auto col = make_col("g", {"A", "B"});
    std::vector<std::string> levels = {};
    REQUIRE_THROWS_AS(encode(col, levels), gelex::InvalidInputException);
}

// ================================================================
// collect_levels
// ================================================================

TEST_CASE("collect_levels returns sorted unique values", "[encode][dataframe]")
{
    auto col = make_col("g", {"B", "A", "B", "C", "A"});
    auto levels = collect_levels(col);

    REQUIRE(levels.size() == 3);
    REQUIRE(levels[0] == "A");
    REQUIRE(levels[1] == "B");
    REQUIRE(levels[2] == "C");
}

TEST_CASE("collect_levels single level does not throw", "[encode][dataframe]")
{
    auto col = make_col("g", {"X", "X", "X"});
    auto levels = collect_levels(col);

    REQUIRE(levels.size() == 1);
    REQUIRE(levels[0] == "X");
}

// ================================================================
// check_levels
// ================================================================

TEST_CASE("check_levels exact match returns ok", "[encode][dataframe]")
{
    auto col = make_col("g", {"A", "B", "C"});
    std::vector<std::string> levels = {"A", "B", "C"};
    auto mismatch = check_levels(col, levels);

    REQUIRE(mismatch.ok());
    REQUIRE(mismatch.missing_in_data.empty());
    REQUIRE(mismatch.missing_in_levels.empty());
}

TEST_CASE(
    "check_levels missing_in_data when levels not in data",
    "[encode][dataframe]")
{
    auto col = make_col("g", {"A", "B"});
    std::vector<std::string> levels = {"A", "B", "Z"};
    auto mismatch = check_levels(col, levels);

    REQUIRE_FALSE(mismatch.ok());
    REQUIRE(mismatch.missing_in_data.size() == 1);
    REQUIRE(mismatch.missing_in_data[0] == "Z");
    REQUIRE(mismatch.missing_in_levels.empty());
}

TEST_CASE(
    "check_levels missing_in_levels when data has unspecified values",
    "[encode][dataframe]")
{
    auto col = make_col("g", {"A", "B", "D"});
    std::vector<std::string> levels = {"A", "B"};
    auto mismatch = check_levels(col, levels);

    REQUIRE_FALSE(mismatch.ok());
    REQUIRE(mismatch.missing_in_levels.size() == 1);
    REQUIRE(mismatch.missing_in_levels[0] == "D");
    REQUIRE(mismatch.missing_in_data.empty());
}

TEST_CASE(
    "check_levels both missing_in_data and missing_in_levels",
    "[encode][dataframe]")
{
    auto col = make_col("g", {"A", "B", "D"});
    std::vector<std::string> levels = {"A", "B", "Z"};
    auto mismatch = check_levels(col, levels);

    REQUIRE_FALSE(mismatch.ok());
    REQUIRE(mismatch.missing_in_data.size() == 1);
    REQUIRE(mismatch.missing_in_data[0] == "Z");
    REQUIRE(mismatch.missing_in_levels.size() == 1);
    REQUIRE(mismatch.missing_in_levels[0] == "D");
}
