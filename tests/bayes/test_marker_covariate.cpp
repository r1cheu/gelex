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
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include "gelex/bayes/marker_covariate.h"
#include "gelex/data/dataframe/column.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/reader.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"

#include "file_fixture.h"

namespace
{
auto read_annotation_frame(
    const std::filesystem::path& path,
    std::size_t annotation_count) -> gelex::DataFrame<std::string>
{
    gelex::ReadOptions options;
    options.delimiter = '\t';
    options.index_cols = {1};
    std::vector<gelex::ColumnType> schema{
        gelex::ColumnType::String,
        gelex::ColumnType::Int,
        gelex::ColumnType::String,
        gelex::ColumnType::String};
    schema.resize(schema.size() + annotation_count, gelex::ColumnType::Double);
    return gelex::read_dataframe<std::string>(path, options, schema);
}
}  // namespace

TEST_CASE(
    "MarkerCovariate aligns marker annotations",
    "[bayes][marker_covariate]")
{
    gelex::test::FileFixture files;
    const auto bim_path
        = files.create_text_file("1 rs1 0 100 A G\n2 rs2 0 200 C T\n", ".bim");
    const auto annotation_path = files.create_text_file(
        "CHR\tSNP\tBP\tA1\tA2\tFrequencyDifference\tCoding\n"
        "2\trs2\t200\tC\tT\t0.75\t1\n"
        "3\textra\t300\tG\tA\t0.50\t0\n"
        "1\trs1\t100\tA\tG\t-0.25\t0\n",
        ".anno");
    const auto marker_metadata = gelex::read_bim(bim_path);
    auto frame = read_annotation_frame(annotation_path, 2);

    const auto marker_covariate = gelex::bayes::make_marker_covariate(
        std::move(frame), marker_metadata);

    const std::vector<std::string> expected_names{
        "Intercept", "FrequencyDifference", "Coding"};
    REQUIRE(
        std::vector(
            marker_covariate.annotation_names().begin(),
            marker_covariate.annotation_names().end())
        == expected_names);
    const Eigen::MatrixXd expected{{1.0, 1.0}, {-0.25, 0.75}, {0.0, 1.0}};
    REQUIRE(marker_covariate.X().isApprox(expected));
}

TEST_CASE(
    "MarkerCovariate requires an annotation column",
    "[bayes][marker_covariate]")
{
    gelex::test::FileFixture files;
    const auto bim_path = files.create_text_file("1 rs1 0 100 A G\n", ".bim");
    const auto annotation_path = files.create_text_file(
        "CHR\tSNP\tBP\tA1\tA2\n1\trs1\t100\tA\tG\n", ".anno");
    const auto marker_metadata = gelex::read_bim(bim_path);
    auto frame = read_annotation_frame(annotation_path, 0);

    REQUIRE_THROWS_AS(
        gelex::bayes::make_marker_covariate(std::move(frame), marker_metadata),
        gelex::GelexException);
}

TEST_CASE(
    "MarkerCovariate requires every target marker",
    "[bayes][marker_covariate]")
{
    gelex::test::FileFixture files;
    const auto bim_path
        = files.create_text_file("1 rs1 0 100 A G\n2 rs2 0 200 C T\n", ".bim");
    const auto annotation_path = files.create_text_file(
        "CHR\tSNP\tBP\tA1\tA2\tAnnotation\n"
        "1\trs1\t100\tA\tG\t0.5\n",
        ".anno");
    const auto marker_metadata = gelex::read_bim(bim_path);
    auto frame = read_annotation_frame(annotation_path, 1);

    REQUIRE_THROWS_AS(
        gelex::bayes::make_marker_covariate(std::move(frame), marker_metadata),
        gelex::GelexException);
}

TEST_CASE(
    "MarkerCovariate validates aligned marker metadata",
    "[bayes][marker_covariate]")
{
    gelex::test::FileFixture files;
    const auto bim_path = files.create_text_file("1 rs1 0 100 A G\n", ".bim");
    const auto marker_metadata = gelex::read_bim(bim_path);

    const auto mismatch_throws = [&](const std::string& row)
    {
        const auto path = files.create_text_file(
            "CHR\tSNP\tBP\tA1\tA2\tAnnotation\n" + row, ".anno");
        auto frame = read_annotation_frame(path, 1);
        REQUIRE_THROWS_AS(
            gelex::bayes::make_marker_covariate(
                std::move(frame), marker_metadata),
            gelex::GelexException);
    };

    mismatch_throws("2\trs1\t100\tA\tG\t0.5\n");
    mismatch_throws("1\trs1\t101\tA\tG\t0.5\n");
    mismatch_throws("1\trs1\t100\tG\tG\t0.5\n");
    mismatch_throws("1\trs1\t100\tA\tT\t0.5\n");
}
