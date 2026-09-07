// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <string>
#include <string_view>

#include "gelex/data/covariates.h"
#include "gelex/data/dataframe/reader.h"

#include "file_fixture.h"

using gelex::read_dataframe;
using gelex::ReadOptions;
using gelex::test::FileFixture;

TEST_CASE(
    "make_quantitative_covariate converts dataframe columns",
    "[data][covariates]")
{
    FileFixture files;
    constexpr std::string_view content
        = "id\theight\tage\n"
          "s1\t1.5\t20\n"
          "s2\t2.5\t30\n";

    auto path = files.create_text_file(content, ".tsv");
    ReadOptions options;
    options.index_cols = {0};
    auto frame = read_dataframe<std::string, double>(path.string(), options);

    auto covariate = gelex::make_quantitative_covariate(frame);

    REQUIRE(covariate.names.size() == 2);
    REQUIRE(covariate.names[0] == "height");
    REQUIRE(covariate.names[1] == "age");

    Eigen::MatrixXd expected{{1.5, 20.0}, {2.5, 30.0}};
    REQUIRE(covariate.X.isApprox(expected));
}

TEST_CASE(
    "make_discrete_covariate dummy encodes multi-level columns",
    "[data][covariates]")
{
    FileFixture files;
    constexpr std::string_view content
        = "id\tsex\tbatch\tmono\n"
          "s1\tF\tB\tz\n"
          "s2\tM\tA\tz\n"
          "s3\tM\tB\tz\n";

    auto path = files.create_text_file(content, ".tsv");
    ReadOptions options;
    options.index_cols = {0};
    auto frame
        = read_dataframe<std::string, std::string>(path.string(), options);

    auto covariate = gelex::make_discrete_covariate(frame);

    REQUIRE(covariate.terms.size() == 2);
    REQUIRE(covariate.terms[0].name == "sex");
    REQUIRE(covariate.terms[1].name == "batch");
    REQUIRE(covariate.terms[0].reference_level == "F");
    REQUIRE(covariate.terms[1].reference_level == "A");

    Eigen::MatrixXd expected{{0.0, 1.0}, {1.0, 0.0}, {1.0, 1.0}};
    REQUIRE(covariate.X.isApprox(expected));
}
