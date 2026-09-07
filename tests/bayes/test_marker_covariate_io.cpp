// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <string>
#include <utility>
#include <vector>

#include "gelex/bayes/marker_covariate.h"
#include "gelex/bayes/marker_covariate_io.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"

#include "file_fixture.h"

TEST_CASE(
    "read_marker_annotation reads the strict annotation format",
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

    auto annotation = gelex::bayes::read_marker_annotation(annotation_path);
    REQUIRE(annotation.rows() == 2);
    REQUIRE(annotation.index().keys()[0] == "rs2");
    REQUIRE(annotation["FrequencyDifference"].as<double>()[0] == 0.75);

    const auto marker_covariate = gelex::bayes::make_marker_covariate(
        std::move(annotation), marker_metadata);

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
    "read_marker_annotation requires the fixed metadata header",
    "[bayes][marker_covariate][io]")
{
    gelex::test::FileFixture files;
    const auto annotation_path = files.create_text_file(
        "CHROM\tSNP\tBP\tA1\tA2\tAnnotation\n"
        "1\trs1\t100\tA\tG\t0.5\n",
        ".anno");

    REQUIRE_THROWS_AS(
        gelex::bayes::read_marker_annotation(annotation_path),
        gelex::GelexException);
}

TEST_CASE(
    "read_marker_annotation rejects duplicate SNPs",
    "[bayes][marker_covariate][io]")
{
    gelex::test::FileFixture files;
    const auto annotation_path = files.create_text_file(
        "CHR\tSNP\tBP\tA1\tA2\tAnnotation\n"
        "1\trs1\t100\tA\tG\t0.5\n"
        "1\trs1\t100\tA\tG\t0.6\n",
        ".anno");

    REQUIRE_THROWS_AS(
        gelex::bayes::read_marker_annotation(annotation_path),
        gelex::GelexException);
}

TEST_CASE(
    "read_marker_annotation rejects non-numeric annotation values",
    "[bayes][marker_covariate][io]")
{
    gelex::test::FileFixture files;
    const auto annotation_path = files.create_text_file(
        "CHR\tSNP\tBP\tA1\tA2\tAnnotation\n"
        "1\trs1\t100\tA\tG\tnot-a-number\n",
        ".anno");

    REQUIRE_THROWS_AS(
        gelex::bayes::read_marker_annotation(annotation_path),
        gelex::GelexException);
}
