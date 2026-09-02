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
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include "gelex/bayes/design.h"
#include "gelex/data/dataframe/constants.h"
#include "gelex/data/dataframe/reader.h"
#include "gelex/exception.h"

#include "file_fixture.h"

static_assert(!std::is_constructible_v<
              gelex::bayes::RandomDesign,
              std::string,
              std::vector<std::string>,
              Eigen::MatrixXd>);

TEST_CASE(
    "make_random_designs builds one-hot coefficient blocks",
    "[bayes][design_factory]")
{
    gelex::test::FileFixture files;
    constexpr std::string_view content
        = "id\tsex\tbatch\n"
          "s1\tF\tB\n"
          "s2\tM\tA\n"
          "s3\tM\tB\n";

    gelex::ReadOptions options;
    options.index_cols = {0};
    const auto frame = gelex::read_dataframe<std::string, std::string>(
        files.create_text_file(content, ".tsv").string(), options);

    const auto designs = gelex::bayes::make_random_designs(frame);

    REQUIRE(designs.size() == 2);
    REQUIRE(designs[0].name() == "sex");
    const std::vector<std::string> expected_sex_names{
        std::string{"sex"} + gelex::separator + "F",
        std::string{"sex"} + gelex::separator + "M"};
    REQUIRE(std::ranges::equal(designs[0].column_names(), expected_sex_names));
    REQUIRE(
        designs[0].X().isApprox(
            Eigen::MatrixXd{{1.0, 0.0}, {0.0, 1.0}, {0.0, 1.0}}));
    REQUIRE(designs[0].xtx_diag().isApprox(Eigen::VectorXd{{1.0, 2.0}}));

    REQUIRE(designs[1].name() == "batch");
    const std::vector<std::string> expected_batch_names{
        std::string{"batch"} + gelex::separator + "A",
        std::string{"batch"} + gelex::separator + "B"};
    REQUIRE(
        std::ranges::equal(designs[1].column_names(), expected_batch_names));
    REQUIRE(
        designs[1].X().isApprox(
            Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {0.0, 1.0}}));
    REQUIRE(designs[1].xtx_diag().isApprox(Eigen::VectorXd{{1.0, 2.0}}));
}

TEST_CASE(
    "make_quantitative_random_design keeps one multi-column block",
    "[bayes][design_factory]")
{
    gelex::test::FileFixture files;
    constexpr std::string_view content
        = "id\tx1\tx2\n"
          "s1\t1.0\t0.0\n"
          "s2\t0.0\t2.0\n"
          "s3\t1.0\t1.0\n";

    gelex::ReadOptions options;
    options.index_cols = {0};
    const auto frame = gelex::read_dataframe<std::string, double>(
        files.create_text_file(content, ".tsv").string(), options);

    const auto design
        = gelex::bayes::make_quantitative_random_design(frame, "qrand");

    REQUIRE(design.name() == "qrand");
    REQUIRE(
        std::ranges::equal(
            design.column_names(), std::vector<std::string>{"x1", "x2"}));
    REQUIRE(design.X().isApprox(
        Eigen::MatrixXd{{1.0, 0.0}, {0.0, 2.0}, {1.0, 1.0}}));
    REQUIRE(design.xtx_diag().isApprox(Eigen::VectorXd{{2.0, 5.0}}));
}

TEST_CASE(
    "RandomDesign factories reject invalid blocks",
    "[bayes][design_factory]")
{
    gelex::test::FileFixture files;

    SECTION("single-level discrete column")
    {
        constexpr std::string_view content
            = "id\tmono\n"
              "s1\tz\n"
              "s2\tz\n";
        gelex::ReadOptions options;
        options.index_cols = {0};
        const auto frame = gelex::read_dataframe<std::string, std::string>(
            files.create_text_file(content, ".tsv").string(), options);

        REQUIRE_THROWS_AS(
            gelex::bayes::make_random_designs(frame), gelex::GelexException);
    }

    SECTION("empty block name")
    {
        constexpr std::string_view content
            = "id\tx\n"
              "s1\t1.0\n"
              "s2\t2.0\n";
        gelex::ReadOptions options;
        options.index_cols = {0};
        const auto frame = gelex::read_dataframe<std::string, double>(
            files.create_text_file(content, ".tsv").string(), options);

        REQUIRE_THROWS_AS(
            gelex::bayes::make_quantitative_random_design(frame, ""),
            gelex::GelexException);
    }

    SECTION("zero quantitative column")
    {
        constexpr std::string_view content
            = "id\tx\n"
              "s1\t0.0\n"
              "s2\t0.0\n";
        gelex::ReadOptions options;
        options.index_cols = {0};
        const auto frame = gelex::read_dataframe<std::string, double>(
            files.create_text_file(content, ".tsv").string(), options);

        REQUIRE_THROWS_AS(
            gelex::bayes::make_quantitative_random_design(frame, "zero"),
            gelex::GelexException);
    }

    SECTION("non-finite column squared norm")
    {
        constexpr std::string_view content
            = "id\tx\n"
              "s1\t1e308\n"
              "s2\t1e308\n";
        gelex::ReadOptions options;
        options.index_cols = {0};
        const auto frame = gelex::read_dataframe<std::string, double>(
            files.create_text_file(content, ".tsv").string(), options);

        REQUIRE_THROWS_AS(
            gelex::bayes::make_quantitative_random_design(frame, "overflow"),
            gelex::GelexException);
    }
}
