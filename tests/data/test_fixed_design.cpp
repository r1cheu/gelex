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

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <fmt/format.h>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

#include "gelex/data/covariates.h"
#include "gelex/data/dataframe/constants.h"
#include "gelex/data/dataframe/reader.h"
#include "gelex/data/fixed_design.h"
#include "gelex/exception.h"

#include "file_fixture.h"

TEST_CASE(
    "FixedDesign names one column per non-reference level",
    "[data][fixed_design]")
{
    gelex::test::FileFixture files;
    constexpr std::string_view qcontent
        = "id\tage\n"
          "s1\t10.0\n"
          "s2\t20.0\n"
          "s3\t15.0\n"
          "s4\t25.0\n"
          "s5\t35.0\n";
    constexpr std::string_view dcontent
        = "id\tsex\tbatch\n"
          "s1\tF\tA\n"
          "s2\tM\tA\n"
          "s3\tF\tB\n"
          "s4\tF\tC\n"
          "s5\tM\tB\n";

    gelex::ReadOptions options;
    options.index_cols = {0};
    auto qframe = gelex::read_dataframe<std::string, double>(
        files.create_text_file(qcontent, ".tsv").string(), options);
    auto dframe = gelex::read_dataframe<std::string, std::string>(
        files.create_text_file(dcontent, ".tsv").string(), options);

    auto design = gelex::FixedDesign::make(
        gelex::make_quantitative_covariate(qframe),
        gelex::make_discrete_covariate(dframe));

    REQUIRE(design.X().cols() == 5);
    const std::vector<std::string> expected_column_names{
        std::string{gelex::intercept_name},
        "age",
        fmt::format("sex{}M", gelex::separator),
        fmt::format("batch{}B", gelex::separator),
        fmt::format("batch{}C", gelex::separator)};
    REQUIRE(std::ranges::equal(design.column_names(), expected_column_names));
    REQUIRE(
        std::ranges::equal(
            design.quantitative_names(), std::vector<std::string>{"age"}));
    REQUIRE(design.discrete_terms().size() == 2);
    REQUIRE(design.discrete_terms()[0].name == "sex");
    REQUIRE(design.discrete_terms()[1].name == "batch");
}

TEST_CASE(
    "FixedDesign rejects matrices without full column rank",
    "[data][fixed_design]")
{
    SECTION("zero column")
    {
        REQUIRE_THROWS_AS(
            gelex::FixedDesign::make(
                gelex::QuantitativeCovariate{
                    .names = {"zero"},
                    .X = Eigen::MatrixXd{{0.0}, {0.0}, {0.0}}}),
            gelex::GelexException);
    }

    SECTION("linearly dependent columns")
    {
        REQUIRE_THROWS_AS(
            gelex::FixedDesign::make(
                gelex::QuantitativeCovariate{
                    .names = {"x", "twice_x"},
                    .X = Eigen::MatrixXd{{1.0, 2.0}, {2.0, 4.0}, {3.0, 6.0}}}),
            gelex::GelexException);
    }
}

TEST_CASE("FixedDesign rejects non-finite matrices", "[data][fixed_design]")
{
    REQUIRE_THROWS_AS(
        gelex::FixedDesign::make(
            gelex::QuantitativeCovariate{
                .names = {"x"},
                .X = Eigen::
                    MatrixXd{{1.0}, {std::numeric_limits<double>::infinity()}}}),
        gelex::GelexException);
}
