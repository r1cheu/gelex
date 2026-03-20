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

#include <cstdint>
#include <format>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "file_fixture.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/dataframe_reader.h"
#include "gelex/exception.h"

using gelex::InvalidInputException;
using gelex::df::ColumnType;
using gelex::df::DataFrame;
using gelex::df::kSeparator;
using gelex::df::NaAction;
using gelex::df::read_dataframe;
using gelex::df::ReadOptions;
using gelex::test::FileFixture;

namespace
{

auto make_basic_df(FileFixture& files) -> DataFrame<std::string>
{
    constexpr std::string_view kContent
        = "id\tx\ty\n"
          "s1\t1.5\t2.5\n"
          "s2\t3.0\t4.0\n"
          "s3\t5.0\t6.0\n";

    auto path = files.create_text_file(kContent, ".tsv");
    ReadOptions opts;
    opts.index_cols = {"id"};
    opts.schema = ColumnType::Double;
    return read_dataframe<std::string>(path.string(), opts);
}

}  // namespace

// ================================================================
// DataFrameReader — basic TSV reading
// ================================================================

TEST_CASE(
    "read_dataframe basic TSV with string index yields correct shape",
    "[dataframe]")
{
    FileFixture files;
    auto df = make_basic_df(files);

    REQUIRE(df.rows() == 3);
    REQUIRE(df.cols() == 2);
}

TEST_CASE("read_dataframe basic TSV yields correct column names", "[dataframe]")
{
    FileFixture files;
    auto df = make_basic_df(files);

    auto names = df.names();
    REQUIRE(names.size() == 2);
    REQUIRE(names[0] == "x");
    REQUIRE(names[1] == "y");
}

TEST_CASE(
    "read_dataframe basic TSV yields correct double column values",
    "[dataframe]")
{
    FileFixture files;
    auto df = make_basic_df(files);

    auto x = df.col("x").as<double>();
    REQUIRE(x[0] == 1.5);
    REQUIRE(x[1] == 3.0);
    REQUIRE(x[2] == 5.0);

    auto y = df.col("y").as<double>();
    REQUIRE(y[0] == 2.5);
    REQUIRE(y[1] == 4.0);
    REQUIRE(y[2] == 6.0);
}

TEST_CASE(
    "read_dataframe basic TSV yields correct row positions",
    "[dataframe]")
{
    FileFixture files;
    auto df = make_basic_df(files);

    REQUIRE(df.row_position("s1") == 0);
    REQUIRE(df.row_position("s2") == 1);
    REQUIRE(df.row_position("s3") == 2);
}

// ================================================================
// DataFrameReader — select_cols
// ================================================================

TEST_CASE(
    "read_dataframe with select_cols picks only requested columns",
    "[dataframe]")
{
    FileFixture files;
    constexpr std::string_view kContent
        = "id\ta\tb\tc\td\n"
          "r1\t1\t2\t3\t4\n"
          "r2\t5\t6\t7\t8\n";

    auto path = files.create_text_file(kContent, ".tsv");
    ReadOptions opts;
    opts.index_cols = {"id"};
    opts.select_cols = {"b", "d"};
    opts.schema = ColumnType::Int;
    auto df = read_dataframe<std::string>(path.string(), opts);

    REQUIRE(df.cols() == 2);
    auto names = df.names();
    REQUIRE(names[0] == "b");
    REQUIRE(names[1] == "d");

    auto b = df.col("b").as<std::int32_t>();
    REQUIRE(b[0] == 2);
    REQUIRE(b[1] == 6);

    auto d = df.col("d").as<std::int32_t>();
    REQUIRE(d[0] == 4);
    REQUIRE(d[1] == 8);
}

// ================================================================
// DataFrameReader — no-header mode
// ================================================================

TEST_CASE(
    "read_dataframe without header generates numeric column names",
    "[dataframe]")
{
    FileFixture files;
    constexpr std::string_view kContent
        = "10\t1.0\n"
          "20\t2.0\n"
          "30\t3.0\n";

    auto path = files.create_text_file(kContent, ".tsv");
    ReadOptions opts;
    opts.header = false;
    opts.schema = std::vector{ColumnType::Int, ColumnType::Double};
    auto df = read_dataframe<std::int32_t>(path.string(), opts);

    REQUIRE(df.rows() == 3);
    REQUIRE(df.cols() == 2);

    auto names = df.names();
    REQUIRE(names[0] == "0");
    REQUIRE(names[1] == "1");
}

// ================================================================
// DataFrameReader — auto-generated integer index
// ================================================================

TEST_CASE(
    "read_dataframe without index_cols generates sequential int index",
    "[dataframe]")
{
    FileFixture files;
    constexpr std::string_view kContent
        = "x\ty\n"
          "1.0\t2.0\n"
          "3.0\t4.0\n"
          "5.0\t6.0\n";

    auto path = files.create_text_file(kContent, ".tsv");
    ReadOptions opts;
    opts.schema = ColumnType::Double;
    auto df = read_dataframe<std::int32_t>(path.string(), opts);

    REQUIRE(df.rows() == 3);
    REQUIRE(df.row_position(0) == 0);
    REQUIRE(df.row_position(1) == 1);
    REQUIRE(df.row_position(2) == 2);
}

// ================================================================
// DataFrameReader — composite index
// ================================================================

TEST_CASE(
    "read_dataframe with two index_cols builds composite string key",
    "[dataframe]")
{
    FileFixture files;
    constexpr std::string_view kContent
        = "fid\tiid\tscore\n"
          "F1\tI1\t0.1\n"
          "F1\tI2\t0.2\n"
          "F2\tI1\t0.3\n";

    auto path = files.create_text_file(kContent, ".tsv");
    ReadOptions opts;
    opts.index_cols = {"fid", "iid"};
    opts.schema = ColumnType::Double;
    auto df = read_dataframe<std::string>(path.string(), opts);

    REQUIRE(df.rows() == 3);

    auto key0 = std::format("{}{}{}", "F1", kSeparator, "I1");
    auto key1 = std::format("{}{}{}", "F1", kSeparator, "I2");
    auto key2 = std::format("{}{}{}", "F2", kSeparator, "I1");

    REQUIRE(df.row_position(key0) == 0);
    REQUIRE(df.row_position(key1) == 1);
    REQUIRE(df.row_position(key2) == 2);
}

// ================================================================
// DataFrameReader — NA handling
// ================================================================

TEST_CASE(
    "read_dataframe with NaAction::Exclude skips rows containing NA",
    "[dataframe]")
{
    FileFixture files;
    constexpr std::string_view kContent
        = "id\tx\n"
          "s1\t1.0\n"
          "s2\tNA\n"
          "s3\t3.0\n";

    auto path = files.create_text_file(kContent, ".tsv");
    ReadOptions opts;
    opts.index_cols = {"id"};
    opts.schema = ColumnType::Double;
    opts.na_action = NaAction::Exclude;
    auto df = read_dataframe<std::string>(path.string(), opts);

    REQUIRE(df.rows() == 2);
    REQUIRE(df.row_position("s1") == 0);
    REQUIRE(df.row_position("s3") == 1);
}

TEST_CASE(
    "read_dataframe with NaAction::Throw raises InvalidInputException on NA",
    "[dataframe]")
{
    FileFixture files;
    constexpr std::string_view kContent
        = "id\tx\n"
          "s1\t1.0\n"
          "s2\tNA\n";

    auto path = files.create_text_file(kContent, ".tsv");
    ReadOptions opts;
    opts.index_cols = {"id"};
    opts.schema = ColumnType::Double;
    opts.na_action = NaAction::Throw;

    REQUIRE_THROWS_AS(
        read_dataframe<std::string>(path.string(), opts),
        InvalidInputException);
}

// ================================================================
// DataFrameReader — parse error
// ================================================================

TEST_CASE(
    "read_dataframe with non-numeric value in double column throws",
    "[dataframe]")
{
    FileFixture files;
    constexpr std::string_view kContent
        = "id\tx\n"
          "s1\tabc\n";

    auto path = files.create_text_file(kContent, ".tsv");
    ReadOptions opts;
    opts.index_cols = {"id"};
    opts.schema = ColumnType::Double;

    REQUIRE_THROWS_AS(
        read_dataframe<std::string>(path.string(), opts),
        InvalidInputException);
}

// ================================================================
// DataFrameReader — comma delimiter
// ================================================================

TEST_CASE(
    "read_dataframe with comma delimiter parses CSV correctly",
    "[dataframe]")
{
    FileFixture files;
    constexpr std::string_view kContent
        = "id,x,y\n"
          "s1,1.5,2.5\n"
          "s2,3.0,4.0\n";

    auto path = files.create_text_file(kContent, ".csv");
    ReadOptions opts;
    opts.delimiter = ',';
    opts.index_cols = {"id"};
    opts.schema = ColumnType::Double;
    auto df = read_dataframe<std::string>(path.string(), opts);

    REQUIRE(df.rows() == 2);
    REQUIRE(df.col("x").as<double>()[0] == 1.5);
    REQUIRE(df.col("y").as<double>()[1] == 4.0);
}

// ================================================================
// DataFrame — col error cases
// ================================================================

TEST_CASE(
    "DataFrame col by non-existent name throws InvalidInputException",
    "[dataframe]")
{
    FileFixture files;
    auto df = make_basic_df(files);

    REQUIRE_THROWS_AS(df.col("missing"), InvalidInputException);
}

TEST_CASE(
    "DataFrame col by out-of-range index throws InvalidInputException",
    "[dataframe]")
{
    FileFixture files;
    auto df = make_basic_df(files);

    REQUIRE_THROWS_AS(df.col(99), InvalidInputException);
}

// ================================================================
// DataFrame — clone
// ================================================================

TEST_CASE("DataFrame clone produces independent copy", "[dataframe]")
{
    FileFixture files;
    auto df = make_basic_df(files);
    auto cloned = df.clone();

    // Mutate clone's column data
    cloned.col("x").as<double>()[0] = 999.0;

    // Original must be unchanged
    REQUIRE(df.col("x").as<double>()[0] == 1.5);
    REQUIRE(cloned.col("x").as<double>()[0] == 999.0);
}

TEST_CASE(
    "DataFrame clone preserves rows, cols, and row positions",
    "[dataframe]")
{
    FileFixture files;
    auto df = make_basic_df(files);
    auto cloned = df.clone();

    REQUIRE(cloned.rows() == df.rows());
    REQUIRE(cloned.cols() == df.cols());
    REQUIRE(cloned.row_position("s2") == df.row_position("s2"));
}

// ================================================================
// DataFrame — gather
// ================================================================

TEST_CASE(
    "DataFrame gather reorders rows and updates row_position",
    "[dataframe]")
{
    FileFixture files;
    auto df = make_basic_df(files);

    // Reorder to {s2, s1}
    std::vector<std::size_t> indices = {1, 0};
    df.gather(indices);

    REQUIRE(df.rows() == 2);
    REQUIRE(df.row_position("s2") == 0);
    REQUIRE(df.row_position("s1") == 1);
}

TEST_CASE(
    "DataFrame gather keeps correct column data after reorder",
    "[dataframe]")
{
    FileFixture files;
    auto df = make_basic_df(files);

    std::vector<std::size_t> indices = {1, 0};
    df.gather(indices);

    auto x = df.col("x").as<double>();
    // s2.x==3.0 now at row 0; s1.x==1.5 at row 1
    REQUIRE(x[0] == 3.0);
    REQUIRE(x[1] == 1.5);
}

// ================================================================
// DataFrame — intersect
// ================================================================

TEST_CASE(
    "DataFrame intersect of two overlapping DataFrames keeps common rows",
    "[dataframe]")
{
    FileFixture files;

    constexpr std::string_view kContent1
        = "id\tv\n"
          "s1\t1.0\n"
          "s2\t2.0\n"
          "s3\t3.0\n";

    constexpr std::string_view kContent2
        = "id\tv\n"
          "s2\t20.0\n"
          "s3\t30.0\n"
          "s4\t40.0\n";

    ReadOptions opts;
    opts.index_cols = {"id"};
    opts.schema = ColumnType::Double;

    auto path1 = files.create_text_file(kContent1, ".tsv");
    auto path2 = files.create_text_file(kContent2, ".tsv");
    auto df1 = read_dataframe<std::string>(path1.string(), opts);
    auto df2 = read_dataframe<std::string>(path2.string(), opts);

    DataFrame<std::string>::intersect({&df1, &df2});

    REQUIRE(df1.rows() == 2);
    REQUIRE(df2.rows() == 2);
}

TEST_CASE(
    "DataFrame intersect aligns row data correctly after intersection",
    "[dataframe]")
{
    FileFixture files;

    constexpr std::string_view kContent1
        = "id\tv\n"
          "s1\t1.0\n"
          "s2\t2.0\n"
          "s3\t3.0\n";

    constexpr std::string_view kContent2
        = "id\tv\n"
          "s2\t20.0\n"
          "s3\t30.0\n"
          "s4\t40.0\n";

    ReadOptions opts;
    opts.index_cols = {"id"};
    opts.schema = ColumnType::Double;

    auto path1 = files.create_text_file(kContent1, ".tsv");
    auto path2 = files.create_text_file(kContent2, ".tsv");
    auto df1 = read_dataframe<std::string>(path1.string(), opts);
    auto df2 = read_dataframe<std::string>(path2.string(), opts);

    DataFrame<std::string>::intersect({&df1, &df2});

    // common keys sorted: {s2, s3}
    REQUIRE(df1.row_position("s2") == 0);
    REQUIRE(df1.row_position("s3") == 1);
    REQUIRE(df2.row_position("s2") == 0);
    REQUIRE(df2.row_position("s3") == 1);

    // df1: s2.v==2.0, s3.v==3.0
    auto v1 = df1.col("v").as<double>();
    REQUIRE(v1[0] == 2.0);
    REQUIRE(v1[1] == 3.0);

    // df2: s2.v==20.0, s3.v==30.0
    auto v2 = df2.col("v").as<double>();
    REQUIRE(v2[0] == 20.0);
    REQUIRE(v2[1] == 30.0);
}

// ================================================================
// DataFrame — to_mat
// ================================================================

TEST_CASE("DataFrame to_mat with no args returns all columns", "[dataframe]")
{
    FileFixture files;
    auto df = make_basic_df(files);

    auto mat = df.to_mat();

    REQUIRE(mat.rows() == 3);
    REQUIRE(mat.cols() == 2);
    Eigen::MatrixXd expected(3, 2);
    expected << 1.5, 2.5, 3.0, 4.0, 5.0, 6.0;
    REQUIRE(mat.isApprox(expected));
}

TEST_CASE("DataFrame to_mat with col_indices selects columns", "[dataframe]")
{
    FileFixture files;
    auto df = make_basic_df(files);

    std::vector<std::size_t> indices = {1};
    auto mat = df.to_mat(indices);

    REQUIRE(mat.rows() == 3);
    REQUIRE(mat.cols() == 1);
    REQUIRE(mat(0, 0) == 2.5);
    REQUIRE(mat(1, 0) == 4.0);
    REQUIRE(mat(2, 0) == 6.0);
}

TEST_CASE("DataFrame to_mat with col_names selects columns", "[dataframe]")
{
    FileFixture files;
    auto df = make_basic_df(files);

    std::vector<std::string_view> names = {"y", "x"};
    auto mat = df.to_mat(std::span<const std::string_view>{names});

    REQUIRE(mat.rows() == 3);
    REQUIRE(mat.cols() == 2);
    // columns in requested order: y first, then x
    REQUIRE(mat(0, 0) == 2.5);
    REQUIRE(mat(0, 1) == 1.5);
}

TEST_CASE("DataFrame to_mat with invalid col_name throws", "[dataframe]")
{
    FileFixture files;
    auto df = make_basic_df(files);

    std::vector<std::string_view> names = {"nonexistent"};
    REQUIRE_THROWS_AS(
        df.to_mat(std::span<const std::string_view>{names}),
        InvalidInputException);
}
