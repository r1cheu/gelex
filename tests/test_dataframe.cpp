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

#include <filesystem>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "file_fixture.h"
#include "gelex/data/frame/dataframe.h"
#include "gelex/data/frame/detail/text_utils.h"
#include "gelex/exception.h"

using gelex::compute_common_index_keys;
using gelex::DataFrame;
using gelex::DataFrameLoadPolicy;
using gelex::FileFormatException;
using gelex::FileOpenException;
using gelex::InconsistentColumnCountException;
using gelex::InvalidOperationException;
using gelex::make_sample_id;
using gelex::MissingValueAction;
using gelex::detail::detect_delimiter;
using gelex::detail::split_line_preserve_empty;
using gelex::test::FileFixture;

namespace
{

template <typename ExceptionType, typename T>
auto require_read_throws(const std::filesystem::path& path) -> void
{
    auto read = [&]() { (void)DataFrame<T>::read(path); };

    REQUIRE_THROWS_AS(read(), ExceptionType);
}

}  // namespace

TEST_CASE(
    "DataFrame loads numeric text with missing policies",
    "[data][dataframe]")
{
    FileFixture files;
    auto path = files.create_named_text_file(
        "frame.csv",
        "FID,IID,value_a,value_b\n"
        "f1,s1,1.5,10.0\n"
        "f2,s2,NA,20.0\n"
        "f3,s3,2.0,NA\n");

    DataFrameLoadPolicy policy;
    policy.missing_value_action = MissingValueAction::UseTypeDefault;

    auto frame = DataFrame<double>::read(path, policy);

    REQUIRE(frame.nrows() == 3);
    REQUIRE(frame.ncols() == 2);
    REQUIRE(frame.column(0).name() == "value_a");
    REQUIRE(frame.column(1).name() == "value_b");
    REQUIRE(frame.index_column().data()[0] == make_sample_id("f1", "s1"));
    REQUIRE(frame.index_column().data()[1] == make_sample_id("f2", "s2"));
    REQUIRE(frame.index_column().data()[2] == make_sample_id("f3", "s3"));
    REQUIRE(frame.column(0).data()[0] == 1.5);
    REQUIRE(frame.column(0).data()[1] == 0.0);
    REQUIRE(frame.column(0).data()[2] == 2.0);
    REQUIRE(frame.column(1).data()[0] == 10.0);
    REQUIRE(frame.column(1).data()[1] == 20.0);
    REQUIRE(frame.column(1).data()[2] == 0.0);
}

TEST_CASE("make_sample_id validates and composes keys", "[data][dataframe]")
{
    std::string key = make_sample_id("fam", "id1");
    REQUIRE(key == std::string("fam") + gelex::kSampleIdSeparator + "id1");
    REQUIRE_THROWS_AS(
        make_sample_id("", "id1"), gelex::ArgumentValidationException);
    REQUIRE_THROWS_AS(
        make_sample_id("fam", ""), gelex::ArgumentValidationException);
}

TEST_CASE(
    "split_line_preserve_empty keeps column positions",
    "[data][dataframe]")
{
    auto tokens = split_line_preserve_empty("a,,c,", ',');

    REQUIRE(tokens.size() == 4);
    REQUIRE(tokens[0] == "a");
    REQUIRE(tokens[1].empty());
    REQUIRE(tokens[2] == "c");
    REQUIRE(tokens[3].empty());
}

TEST_CASE(
    "detect_delimiter picks tab then comma then space",
    "[data][dataframe]")
{
    REQUIRE(detect_delimiter("FID\tIID\tV") == '\t');
    REQUIRE(detect_delimiter("FID,IID,V") == ',');
    REQUIRE(detect_delimiter("FID IID V") == ' ');
}

TEST_CASE(
    "DataFrame default policy skips rows with missing values",
    "[data][dataframe]")
{
    FileFixture files;
    auto path = files.create_named_text_file(
        "frame_skip.csv",
        "FID,IID,value_a,value_b\n"
        "f1,s1,1.5,10.0\n"
        "f2,s2,NA,20.0\n"
        "f3,s3,2.0,NA\n"
        "f4,s4,3.0,40.0\n");

    auto frame = DataFrame<double>::read(path);

    REQUIRE(frame.nrows() == 2);
    REQUIRE(frame.ncols() == 2);
    REQUIRE(frame.index_column().data()[0] == make_sample_id("f1", "s1"));
    REQUIRE(frame.index_column().data()[1] == make_sample_id("f4", "s4"));
    REQUIRE(frame.column(0).data()[0] == 1.5);
    REQUIRE(frame.column(0).data()[1] == 3.0);
    REQUIRE(frame.column(1).data()[0] == 10.0);
    REQUIRE(frame.column(1).data()[1] == 40.0);
}

TEST_CASE("DataFrame rejects inconsistent column counts", "[data][dataframe]")
{
    FileFixture files;
    auto path = files.create_named_text_file(
        "bad.csv",
        "FID,IID,value\n"
        "f1,i1,1\n"
        "f2,i2,2,3\n");

    require_read_throws<InconsistentColumnCountException, int>(path);
}

TEST_CASE("DataFrame rejects duplicated index key on load", "[data][dataframe]")
{
    FileFixture files;
    auto path = files.create_named_text_file(
        "dup.csv",
        "FID,IID,value\n"
        "f1,i1,1\n"
        "f1,i1,2\n");

    require_read_throws<InvalidOperationException, int>(path);
}

TEST_CASE(
    "DataFrame requires FID IID and at least one value column",
    "[data][dataframe]")
{
    FileFixture files;
    auto path = files.create_named_text_file(
        "too_few_columns.csv",
        "FID,IID\n"
        "f1,i1\n"
        "f2,i2\n");

    require_read_throws<InconsistentColumnCountException, int>(path);
}

TEST_CASE(
    "DataFrame requires header to start with FID and IID",
    "[data][dataframe]")
{
    FileFixture files;
    auto path = files.create_named_text_file(
        "bad_header_prefix.csv",
        "sample_id,individual_id,value\n"
        "f1,i1,1\n"
        "f2,i2,2\n");

    require_read_throws<FileFormatException, int>(path);
}

TEST_CASE("DataFrame rejects rows with empty FID or IID", "[data][dataframe]")
{
    FileFixture files;

    auto empty_fid = files.create_named_text_file(
        "empty_fid.csv",
        "FID,IID,value\n"
        ",i1,1\n");
    require_read_throws<gelex::DataParseException, int>(empty_fid);

    auto empty_iid = files.create_named_text_file(
        "empty_iid.csv",
        "FID,IID,value\n"
        "f1,,1\n");
    require_read_throws<gelex::DataParseException, int>(empty_iid);
}

TEST_CASE(
    "DataFrame missing policy Throw rejects missing values",
    "[data][dataframe]")
{
    FileFixture files;
    auto path = files.create_named_text_file(
        "missing_throw.csv",
        "FID,IID,value\n"
        "f1,i1,NA\n");

    DataFrameLoadPolicy policy;
    policy.missing_value_action = MissingValueAction::Throw;

    REQUIRE_THROWS_AS(
        DataFrame<double>::read(path, policy), gelex::DataParseException);
}

TEST_CASE("DataFrame rejects invalid numeric tokens", "[data][dataframe]")
{
    FileFixture files;
    auto path = files.create_named_text_file(
        "invalid_number.csv",
        "FID,IID,value\n"
        "f1,i1,abc\n");

    REQUIRE_THROWS_AS(DataFrame<double>::read(path), gelex::DataParseException);
}

TEST_CASE(
    "DataFrame rejects numbers with trailing characters",
    "[data][dataframe]")
{
    FileFixture files;
    auto path = files.create_named_text_file(
        "invalid_number_trailing.csv",
        "FID,IID,value\n"
        "f1,i1,1x\n");

    REQUIRE_THROWS_AS(DataFrame<double>::read(path), gelex::DataParseException);
}

TEST_CASE("DataFrame read throws on missing file", "[data][dataframe]")
{
    std::filesystem::path missing
        = std::filesystem::temp_directory_path() / "gelex_missing_file.csv";

    REQUIRE_THROWS_AS(DataFrame<int>::read(missing), FileOpenException);
}

TEST_CASE("DataFrame read throws on empty file", "[data][dataframe]")
{
    FileFixture files;
    auto path = files.create_named_text_file("empty.csv", "");

    REQUIRE_THROWS_AS(DataFrame<int>::read(path), FileFormatException);
}

TEST_CASE(
    "DataFrame intersect_index_inplace follows key order",
    "[data][dataframe]")
{
    FileFixture files;
    auto path = files.create_named_text_file(
        "stable.csv",
        "FID,IID,value\n"
        "fa,a,1\n"
        "fb,b,2\n"
        "fc,c,3\n"
        "fd,d,4\n");

    auto frame = DataFrame<int>::read(path);

    std::vector<std::string> keys
        = {make_sample_id("fd", "d"),
           make_sample_id("fb", "b"),
           make_sample_id("fx", "x")};
    frame.intersect_index_inplace(keys);

    REQUIRE(frame.nrows() == 2);
    REQUIRE(frame.index_column().data()[0] == make_sample_id("fd", "d"));
    REQUIRE(frame.index_column().data()[1] == make_sample_id("fb", "b"));
    REQUIRE(frame.column(0).data()[0] == 4);
    REQUIRE(frame.column(0).data()[1] == 2);
}

TEST_CASE(
    "DataFrame intersect_index_inplace rejects duplicate keys",
    "[data][dataframe]")
{
    FileFixture files;
    auto path = files.create_named_text_file(
        "dup_keys.csv",
        "FID,IID,value\n"
        "fa,a,1\n"
        "fb,b,2\n");

    auto frame = DataFrame<int>::read(path);

    std::string key = make_sample_id("fa", "a");
    std::vector<std::string> keys = {key, key};
    REQUIRE_THROWS_AS(
        frame.intersect_index_inplace(keys), InvalidOperationException);
}

TEST_CASE("compute_common_index_keys finds shared ids", "[data][dataframe]")
{
    FileFixture files;
    auto path_a = files.create_named_text_file(
        "common_a.csv",
        "FID,IID,value\n"
        "f1,i1,10\n"
        "f2,i2,20\n"
        "f3,i3,30\n");
    auto path_b = files.create_named_text_file(
        "common_b.csv",
        "FID,IID,value\n"
        "f0,i0,1\n"
        "f2,i2,2\n"
        "f3,i3,3\n");
    auto path_c = files.create_named_text_file(
        "common_c.csv",
        "FID,IID,value\n"
        "f3,i3,300\n"
        "f2,i2,200\n"
        "f9,i9,900\n");

    auto frame_a = DataFrame<int>::read(path_a);
    auto frame_b = DataFrame<int>::read(path_b);
    auto frame_c = DataFrame<int>::read(path_c);

    std::vector<const DataFrame<int>*> frames = {&frame_a, &frame_b, &frame_c};
    auto common = compute_common_index_keys<int>(frames);

    REQUIRE(common.size() == 2);
    REQUIRE(common[0] == make_sample_id("f2", "i2"));
    REQUIRE(common[1] == make_sample_id("f3", "i3"));
}

TEST_CASE(
    "compute_common_index_keys handles empty and null pointers",
    "[data][dataframe]")
{
    std::vector<const DataFrame<int>*> empty_frames;
    REQUIRE(compute_common_index_keys<int>(empty_frames).empty());

    std::vector<const DataFrame<int>*> null_first = {nullptr};
    REQUIRE_THROWS_AS(
        compute_common_index_keys<int>(null_first),
        gelex::ArgumentValidationException);

    FileFixture files;
    auto path = files.create_named_text_file(
        "one_frame.csv",
        "FID,IID,value\n"
        "f1,i1,1\n");
    auto frame = DataFrame<int>::read(path);

    std::vector<const DataFrame<int>*> null_later = {&frame, nullptr};
    REQUIRE_THROWS_AS(
        compute_common_index_keys<int>(null_later),
        gelex::ArgumentValidationException);
}

TEST_CASE("DataFrame column access validates bounds", "[data][dataframe]")
{
    FileFixture files;
    auto path = files.create_named_text_file(
        "column_bounds.csv",
        "FID,IID,value\n"
        "f1,i1,1\n");
    auto frame = DataFrame<int>::read(path);

    REQUIRE_THROWS_AS(frame.column(1), gelex::ColumnRangeException);
}

TEST_CASE("DataFrame columns returns data column names", "[data][dataframe]")
{
    FileFixture files;
    auto path = files.create_named_text_file(
        "columns.csv",
        "FID,IID,age,height\n"
        "f1,i1,10,1\n"
        "f2,i2,20,2\n");

    auto frame = DataFrame<int>::read(path);
    auto names = frame.columns();

    REQUIRE(names.size() == 2);
    REQUIRE(names[0] == "age");
    REQUIRE(names[1] == "height");
}

TEST_CASE("DataFrame eigen materializes data matrix", "[data][dataframe]")
{
    FileFixture files;
    auto path = files.create_named_text_file(
        "eigen.csv",
        "FID,IID,age,height\n"
        "f1,i1,10,1.5\n"
        "f2,i2,20,2.5\n"
        "f3,i3,30,3.5\n");

    auto frame = DataFrame<double>::read(path);
    auto matrix = frame.eigen();

    REQUIRE(matrix.rows() == 3);
    REQUIRE(matrix.cols() == 2);
    REQUIRE(matrix(0, 0) == 10.0);
    REQUIRE(matrix(1, 0) == 20.0);
    REQUIRE(matrix(2, 0) == 30.0);
    REQUIRE(matrix(0, 1) == 1.5);
    REQUIRE(matrix(1, 1) == 2.5);
    REQUIRE(matrix(2, 1) == 3.5);
}

TEST_CASE("DataFrame eigen follows intersected row order", "[data][dataframe]")
{
    FileFixture files;
    auto path = files.create_named_text_file(
        "eigen_intersect.csv",
        "FID,IID,value\n"
        "fa,a,1\n"
        "fb,b,2\n"
        "fc,c,3\n");

    auto frame = DataFrame<int>::read(path);
    std::vector<std::string> keys
        = {make_sample_id("fc", "c"), make_sample_id("fa", "a")};
    frame.intersect_index_inplace(keys);

    auto matrix = frame.eigen();
    REQUIRE(matrix.rows() == 2);
    REQUIRE(matrix.cols() == 1);
    REQUIRE(matrix(0, 0) == 3);
    REQUIRE(matrix(1, 0) == 1);
}

TEST_CASE("Column gather_inplace validates row bounds", "[data][dataframe]")
{
    gelex::Column<int> column("value");
    column.append(10);

    std::vector<size_t> keep_rows = {1};
    REQUIRE_THROWS_AS(
        column.gather_inplace(keep_rows), gelex::ColumnRangeException);
}
