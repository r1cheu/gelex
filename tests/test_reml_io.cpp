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
#include <fmt/format.h>
#include <fstream>
#include <iterator>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "gelex/data/dataframe/constants.h"
#include "gelex/exception.h"
#include "gelex/freq/design.h"
#include "gelex/freq/model.h"
#include "gelex/io/reml.h"
#include "gelex/types/fixed_designs.h"

#include "file_fixture.h"
#include "sample_id_fixture.h"

TEST_CASE("REML writers write result files", "[reml][io]")
{
    std::vector<gelex::freq::RandomDesign> random;
    gelex::freq::RandomDesign animal;
    animal.name = "animal";
    animal.K = Eigen::MatrixXd::Identity(2, 2);
    random.push_back(std::move(animal));

    gelex::FreqModel model{
        Eigen::VectorXd{{1.0, 2.0}},
        gelex::FixedDesign::make(2),
        std::move(random)};
    gelex::FreqState state(model);
    state.fixed().coeffs = Eigen::VectorXd{{1.5}};
    state.fixed().se = Eigen::VectorXd{{0.25}};
    state.random()[0].blup = Eigen::VectorXd{{0.25, -0.5}};
    state.random()[0].variance = 2.0;
    state.random()[0].variance_se = 0.3;
    state.random()[0].variance_ratio = 0.4;
    state.random()[0].variance_ratio_se = 0.05;
    state.residual().variance = 3.0;
    state.residual().variance_se = 0.4;

    const std::vector<std::string> sample_ids{
        gelex::make_sample_id("F1", "I1"), gelex::make_sample_id("F2", "I2")};

    SECTION("write_summary writes fixed and variance components")
    {
        gelex::test::FileFixture files;
        const auto prefix = files.get_test_dir() / "reml_result";
        gelex::write_summary(model, state, -123.456, prefix.string());

        auto summary_path = prefix;
        summary_path += ".summary";
        std::ifstream input(summary_path);
        const std::string content{
            std::istreambuf_iterator<char>{input},
            std::istreambuf_iterator<char>{}};

        REQUIRE(
            content.find("term\ttype\testimate\tse\tratio\tratio_se\tpvalue\n")
            == 0);
        REQUIRE(
            content.find(
                "Intercept\tfixed\t1.50000000e+00\t2.50000000e-01\t-\t-")
            != std::string::npos);
        REQUIRE(
            content.find(
                "animal\tvariance\t2.00000000e+00\t3.00000000e-01\t"
                "4.00000000e-01\t5.00000000e-02")
            != std::string::npos);
        REQUIRE(
            content.find(
                "Residual\tvariance\t3.00000000e+00\t4.00000000e-01\t-\t-")
            != std::string::npos);
        REQUIRE(
            content.find("logL\tmodelfit\t-1.23456000e+02\t-\t-\t-")
            != std::string::npos);
    }

    SECTION("write_effects writes fixed, random, and total effects")
    {
        gelex::test::FileFixture files;
        const auto prefix = files.get_test_dir() / "reml_effects";
        gelex::write_effects(model, state, sample_ids, prefix.string());

        auto effects_path = prefix;
        effects_path += ".effects";
        std::ifstream input(effects_path);
        const std::string content{
            std::istreambuf_iterator<char>{input},
            std::istreambuf_iterator<char>{}};

        REQUIRE(content.find("FID\tIID\tIntercept\tanimal\tTOTAL\n") == 0);
        REQUIRE(
            content.find(
                "F1\tI1\t1.50000000e+00\t2.50000000e-01\t1.75000000e+00\n")
            != std::string::npos);
        REQUIRE(
            content.find(
                "F2\tI2\t1.50000000e+00\t-5.00000000e-01\t1.00000000e+00\n")
            != std::string::npos);
    }

    SECTION("write_effects rejects sample ID count mismatch")
    {
        const std::vector<std::string> mismatch{
            gelex::make_sample_id("F1", "I1")};
        gelex::test::FileFixture files;
        const auto prefix = files.get_test_dir() / "reml_bad_effects";

        REQUIRE_THROWS_AS(
            gelex::write_effects(model, state, mismatch, prefix.string()),
            gelex::GelexException);
    }

    SECTION("write_effects handles models without random effects")
    {
        gelex::FreqModel fixed_only_model{
            Eigen::VectorXd{{1.0, 2.0}}, gelex::FixedDesign::make(2), {}};
        gelex::FreqState fixed_only_state(fixed_only_model);
        fixed_only_state.fixed().coeffs = Eigen::VectorXd{{1.25}};
        fixed_only_state.fixed().se = Eigen::VectorXd{{0.2}};

        gelex::test::FileFixture files;
        const auto prefix = files.get_test_dir() / "reml_fixed_only";
        gelex::write_effects(
            fixed_only_model, fixed_only_state, sample_ids, prefix.string());

        auto effects_path = prefix;
        effects_path += ".effects";
        std::ifstream input(effects_path);
        const std::string content{
            std::istreambuf_iterator<char>{input},
            std::istreambuf_iterator<char>{}};

        REQUIRE(content.find("FID\tIID\tIntercept\tTOTAL\n") == 0);
        REQUIRE(
            content.find("F1\tI1\t1.25000000e+00\t1.25000000e+00\n")
            != std::string::npos);
    }

    SECTION("write_effects writes fixed effect contribution columns")
    {
        gelex::FixedDesign fixed{
            .names = {"Intercept", "age", "batch"},
            .levels
            = {std::nullopt, std::nullopt, std::vector<std::string>{"A", "B"}},
            .reference_levels
            = {std::nullopt, std::nullopt, std::string{"base"}},
            .column_names
            = {"Intercept",
               "age",
               fmt::format("batch{}A", gelex::separator),
               fmt::format("batch{}B", gelex::separator)},
            .X = Eigen::MatrixXd{{1.0, 2.0, 1.0, 0.0}, {1.0, 4.0, 0.0, 1.0}},
            .XtX_diag = Eigen::VectorXd{{2.0, 20.0, 1.0, 1.0}}};
        gelex::FreqModel fixed_model{
            Eigen::VectorXd{{1.0, 2.0}}, std::move(fixed), {}};
        gelex::FreqState fixed_state(fixed_model);
        fixed_state.fixed().coeffs = Eigen::VectorXd{{1.0, 0.25, 2.0, 0.75}};
        fixed_state.fixed().se = Eigen::VectorXd{{0.1, 0.2, 0.3, 0.4}};

        gelex::test::FileFixture files;
        const auto prefix = files.get_test_dir() / "reml_fixed_columns";
        gelex::write_effects(
            fixed_model, fixed_state, sample_ids, prefix.string());

        auto effects_path = prefix;
        effects_path += ".effects";
        std::ifstream input(effects_path);
        const std::string content{
            std::istreambuf_iterator<char>{input},
            std::istreambuf_iterator<char>{}};

        REQUIRE(
            content.find(
                fmt::format(
                    "FID\tIID\tIntercept\tage\tbatch{0}A\tbatch{0}B\tTOTAL\n",
                    gelex::separator))
            == 0);
        REQUIRE(
            content.find(
                "F1\tI1\t1.00000000e+00\t5.00000000e-01\t"
                "2.00000000e+00\t0.00000000e+00\t3.50000000e+00\n")
            != std::string::npos);
        REQUIRE(
            content.find(
                "F2\tI2\t1.00000000e+00\t1.00000000e+00\t"
                "0.00000000e+00\t7.50000000e-01\t2.75000000e+00\n")
            != std::string::npos);
    }
}
