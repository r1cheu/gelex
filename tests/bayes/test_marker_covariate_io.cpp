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
#include <string>
#include <vector>

#include "gelex/bayes/marker_covariate_io.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"

#include "file_fixture.h"

TEST_CASE(
    "read_marker_covariate reads the strict annotation format",
    "[bayes][marker_covariate][io]")
{
    gelex::test::FileFixture files;
    const auto bim_path
        = files.create_text_file("1 rs1 0 100 A G\n2 rs2 0 200 C T\n", ".bim");
    const auto annotation_path = files.create_text_file(
        "CHR\tSNP\tBP\tA1\tA2\tFrequencyDifference\n"
        "2\trs2\t200\tC\tT\t0.75\n"
        "1\trs1\t100\tA\tG\t-0.25\n",
        ".anno");
    const auto marker_metadata = gelex::read_bim(bim_path);

    const auto marker_covariate
        = gelex::bayes::read_marker_covariate(annotation_path, marker_metadata);

    const std::vector<std::string> expected_names{
        "Intercept", "FrequencyDifference"};
    REQUIRE(
        std::vector(
            marker_covariate.annotation_names().begin(),
            marker_covariate.annotation_names().end())
        == expected_names);
    const Eigen::MatrixXd expected{{1.0, 1.0}, {-0.25, 0.75}};
    REQUIRE(marker_covariate.X().isApprox(expected));
}

TEST_CASE(
    "read_marker_covariate requires the fixed metadata header",
    "[bayes][marker_covariate][io]")
{
    gelex::test::FileFixture files;
    const auto bim_path = files.create_text_file("1 rs1 0 100 A G\n", ".bim");
    const auto annotation_path = files.create_text_file(
        "CHR\tID\tBP\tA1\tA2\tAnnotation\n"
        "1\trs1\t100\tA\tG\t0.5\n",
        ".anno");
    const auto marker_metadata = gelex::read_bim(bim_path);

    REQUIRE_THROWS_AS(
        gelex::bayes::read_marker_covariate(annotation_path, marker_metadata),
        gelex::GelexException);
}

TEST_CASE(
    "read_marker_covariate rejects duplicate SNPs",
    "[bayes][marker_covariate][io]")
{
    gelex::test::FileFixture files;
    const auto bim_path = files.create_text_file("1 rs1 0 100 A G\n", ".bim");
    const auto annotation_path = files.create_text_file(
        "CHR\tSNP\tBP\tA1\tA2\tAnnotation\n"
        "1\trs1\t100\tA\tG\t0.5\n"
        "1\trs1\t100\tA\tG\t0.6\n",
        ".anno");
    const auto marker_metadata = gelex::read_bim(bim_path);

    REQUIRE_THROWS_AS(
        gelex::bayes::read_marker_covariate(annotation_path, marker_metadata),
        gelex::GelexException);
}

TEST_CASE(
    "read_marker_covariate rejects non-numeric annotation values",
    "[bayes][marker_covariate][io]")
{
    gelex::test::FileFixture files;
    const auto bim_path = files.create_text_file("1 rs1 0 100 A G\n", ".bim");
    const auto annotation_path = files.create_text_file(
        "CHR\tSNP\tBP\tA1\tA2\tAnnotation\n"
        "1\trs1\t100\tA\tG\tnot-a-number\n",
        ".anno");
    const auto marker_metadata = gelex::read_bim(bim_path);

    REQUIRE_THROWS_AS(
        gelex::bayes::read_marker_covariate(annotation_path, marker_metadata),
        gelex::GelexException);
}
