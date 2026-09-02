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
#include <catch2/catch_test_macros.hpp>
#include <string>
#include <string_view>

#include "gelex/data/dataframe/constants.h"
#include "gelex/data/dataframe/reader.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"
#include "gelex/freq/design_factory.h"

#include "file_fixture.h"

TEST_CASE(
    "make_random_designs builds sample relationship matrices",
    "[freq][design_factory]")
{
    gelex::test::FileFixture files;
    constexpr std::string_view content
        = "id\tsex\tbatch\n"
          "s1\tF\tB\n"
          "s2\tM\tA\n"
          "s3\tM\tB\n";

    auto path = files.create_text_file(content, ".tsv");
    gelex::ReadOptions options;
    options.index_cols = {0};
    auto frame = gelex::read_dataframe<std::string, std::string>(
        path.string(), options);

    auto designs = gelex::freq::make_random_designs(frame);

    REQUIRE(designs.size() == 2);
    REQUIRE(designs[0].name == "sex");
    REQUIRE(designs[0].levels.has_value());
    const auto& sex_levels = *designs[0].levels;
    REQUIRE(sex_levels.size() == 2);
    REQUIRE(sex_levels[0] == std::string{"sex"} + gelex::separator + "F");
    REQUIRE(sex_levels[1] == std::string{"sex"} + gelex::separator + "M");

    Eigen::MatrixXd expected_sex{
        {1.0, 0.0, 0.0}, {0.0, 1.0, 1.0}, {0.0, 1.0, 1.0}};
    REQUIRE(designs[0].K.isApprox(expected_sex));

    REQUIRE(designs[1].name == "batch");
    REQUIRE(designs[1].levels.has_value());
    const auto& batch_levels = *designs[1].levels;
    REQUIRE(batch_levels.size() == 2);
    REQUIRE(batch_levels[0] == std::string{"batch"} + gelex::separator + "A");
    REQUIRE(batch_levels[1] == std::string{"batch"} + gelex::separator + "B");

    Eigen::MatrixXd expected_batch{
        {1.0, 0.0, 1.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 1.0}};
    REQUIRE(designs[1].K.isApprox(expected_batch));
}

TEST_CASE(
    "make_quantitative_random_design builds linear kernel ZZ^T",
    "[freq][design_factory]")
{
    gelex::test::FileFixture files;
    constexpr std::string_view content
        = "FID\tIID\tx1\tx2\n"
          "F1\ts1\t1.0\t0.0\n"
          "F1\ts2\t0.0\t2.0\n"
          "F1\ts3\t1.0\t1.0\n";

    auto path = files.create_text_file(content, ".tsv");
    auto frame = gelex::read_qcovar(path);

    auto design = gelex::freq::make_quantitative_random_design(frame, "qrand");

    REQUIRE(design.name == "qrand");
    REQUIRE_FALSE(design.levels.has_value());
    REQUIRE_FALSE(design.Z.has_value());

    Eigen::MatrixXd expected{{1.0, 0.0, 1.0}, {0.0, 4.0, 2.0}, {1.0, 2.0, 2.0}};
    REQUIRE(design.K.isApprox(expected));
}

TEST_CASE(
    "make_interaction_design forms rescaled Hadamard product",
    "[freq][design_factory]")
{
    Eigen::MatrixXd lhs{{2.0, 1.0}, {1.0, 4.0}};
    Eigen::MatrixXd rhs{{1.0, 0.5}, {0.5, 3.0}};

    auto design = gelex::freq::make_interaction_design("A:B", lhs, rhs);

    REQUIRE(design.name == "A:B");
    REQUIRE(design.kind == gelex::freq::RandomKind::Interaction);
    REQUIRE_FALSE(design.Z.has_value());

    Eigen::MatrixXd hadamard = lhs.cwiseProduct(rhs);
    Eigen::MatrixXd expected = hadamard / hadamard.diagonal().mean();
    REQUIRE(design.K.isApprox(expected));
}

TEST_CASE(
    "make_random_designs rejects single-level columns",
    "[freq][design_factory]")
{
    gelex::test::FileFixture files;
    constexpr std::string_view content
        = "id\tmono\n"
          "s1\tz\n"
          "s2\tz\n";

    auto path = files.create_text_file(content, ".tsv");
    gelex::ReadOptions options;
    options.index_cols = {0};
    auto frame = gelex::read_dataframe<std::string, std::string>(
        path.string(), options);

    REQUIRE_THROWS_AS(
        gelex::freq::make_random_designs(frame), gelex::GelexException);
}
