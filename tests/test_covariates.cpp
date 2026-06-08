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

#include <string>
#include <string_view>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "file_fixture.h"
#include "gelex/data/covariates.h"
#include "gelex/data/dataframe/constants.h"
#include "gelex/data/dataframe/reader.h"
#include "gelex/exception.h"

using gelex::dataframe::read_dataframe;
using gelex::dataframe::ReadOptions;
using gelex::dataframe::SEPARATOR;
using gelex::test::FileFixture;

TEST_CASE(
    "make_quantitative_covariate converts dataframe columns",
    "[covariates]")
{
    FileFixture files;
    constexpr std::string_view CONTENT
        = "id\theight\tage\n"
          "s1\t1.5\t20\n"
          "s2\t2.5\t30\n";

    auto path = files.create_text_file(CONTENT, ".tsv");
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
    "[covariates]")
{
    FileFixture files;
    constexpr std::string_view CONTENT
        = "id\tsex\tbatch\tmono\n"
          "s1\tF\tB\tz\n"
          "s2\tM\tA\tz\n"
          "s3\tM\tB\tz\n";

    auto path = files.create_text_file(CONTENT, ".tsv");
    ReadOptions options;
    options.index_cols = {0};
    auto frame
        = read_dataframe<std::string, std::string>(path.string(), options);

    auto covariate = gelex::make_discrete_covariate(frame);

    REQUIRE(covariate.names.size() == 2);
    REQUIRE(covariate.names[0] == "sex");
    REQUIRE(covariate.names[1] == "batch");
    REQUIRE(covariate.reference_levels[0] == "F");
    REQUIRE(covariate.reference_levels[1] == "A");

    Eigen::MatrixXd expected{{0.0, 1.0}, {1.0, 0.0}, {1.0, 1.0}};
    REQUIRE(covariate.X.isApprox(expected));
}

TEST_CASE(
    "make_random_designs builds sample relationship matrices",
    "[covariates]")
{
    FileFixture files;
    constexpr std::string_view CONTENT
        = "id\tsex\tbatch\n"
          "s1\tF\tB\n"
          "s2\tM\tA\n"
          "s3\tM\tB\n";

    auto path = files.create_text_file(CONTENT, ".tsv");
    ReadOptions options;
    options.index_cols = {0};
    auto frame
        = read_dataframe<std::string, std::string>(path.string(), options);

    auto designs = gelex::make_random_designs(frame);

    REQUIRE(designs.size() == 2);
    REQUIRE(designs[0].name == "sex");
    REQUIRE(designs[0].levels.size() == 2);
    REQUIRE(designs[0].levels[0] == std::string{"sex"} + SEPARATOR + "F");
    REQUIRE(designs[0].levels[1] == std::string{"sex"} + SEPARATOR + "M");

    Eigen::MatrixXd expected_sex{
        {1.0, 0.0, 0.0}, {0.0, 1.0, 1.0}, {0.0, 1.0, 1.0}};
    REQUIRE(designs[0].K.isApprox(expected_sex));

    REQUIRE(designs[1].name == "batch");
    REQUIRE(designs[1].levels.size() == 2);
    REQUIRE(designs[1].levels[0] == std::string{"batch"} + SEPARATOR + "A");
    REQUIRE(designs[1].levels[1] == std::string{"batch"} + SEPARATOR + "B");

    Eigen::MatrixXd expected_batch{
        {1.0, 0.0, 1.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 1.0}};
    REQUIRE(designs[1].K.isApprox(expected_batch));
}

TEST_CASE("make_random_designs rejects single-level columns", "[covariates]")
{
    FileFixture files;
    constexpr std::string_view CONTENT
        = "id\tmono\n"
          "s1\tz\n"
          "s2\tz\n";

    auto path = files.create_text_file(CONTENT, ".tsv");
    ReadOptions options;
    options.index_cols = {0};
    auto frame
        = read_dataframe<std::string, std::string>(path.string(), options);

    REQUIRE_THROWS_AS(gelex::make_random_designs(frame), gelex::GelexException);
}
