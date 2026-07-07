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

#include <optional>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <string_view>

#include "cli/predict/compute.h"
#include "file_fixture.h"
#include "gelex/data/dataframe/constants.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/reader.h"

using gelex::SEPARATOR;
using gelex::test::FileFixture;

TEST_CASE(
    "read_param reads term names and coefficients from TSV",
    "[predict][io]")
{
    FileFixture files;
    constexpr std::string_view CONTENT
        = "term\tmean\n"
          "Intercept\t0.5\n"
          "Age\t-0.3\n"
          "Sex_M\t1.2\n";

    auto path = files.create_text_file(CONTENT, ".tsv");
    auto param = gelex::read_param(path);

    auto names = param.index().keys();
    REQUIRE(names.size() == 3);
    REQUIRE(names[0] == "Intercept");
    REQUIRE(names[1] == "Age");
    REQUIRE(names[2] == "Sex_M");

    Eigen::VectorXd values = param["mean"].to_map<double>();
    REQUIRE(values.size() == 3);
    REQUIRE(values[0] == 0.5);
    REQUIRE(values[1] == -0.3);
    REQUIRE(values[2] == 1.2);
}

TEST_CASE(
    "build_covariate_design builds design matrix from qcovar and dcovar",
    "[predict][io]")
{
    FileFixture files;
    auto sep = std::string(1, SEPARATOR);

    // qcovar file (tab-delimited, header)
    auto qcovar_path = files.create_text_file(
        "FID\tIID\tAge\n"
        "F1\tI1\t25.0\n"
        "F1\tI2\t30.0\n"
        "F2\tI1\t35.0\n",
        ".tsv");

    // dcovar file (tab-delimited, header)
    auto dcovar_path = files.create_text_file(
        "FID\tIID\tSex\n"
        "F1\tI1\tM\n"
        "F1\tI2\tF\n"
        "F2\tI1\tM\n",
        ".tsv");

    std::optional<gelex::DataFrame<std::string>> qcovar_df
        = gelex::read_qcovar(qcovar_path);
    std::optional<gelex::DataFrame<std::string>> dcovar_df
        = gelex::read_dcovar(dcovar_path);

    std::vector<std::string> term_names{"Intercept", "Age", "Sex" + sep + "M"};

    auto covariates
        = cli::build_covariate_design(term_names, qcovar_df, dcovar_df, 3);

    REQUIRE(covariates.rows() == 3);
    REQUIRE(covariates.cols() == 3);

    // col 0: Intercept (all 1s)
    REQUIRE(covariates(0, 0) == 1.0);
    REQUIRE(covariates(1, 0) == 1.0);
    REQUIRE(covariates(2, 0) == 1.0);

    // col 1: Age
    REQUIRE(covariates(0, 1) == 25.0);
    REQUIRE(covariates(1, 1) == 30.0);
    REQUIRE(covariates(2, 1) == 35.0);

    // col 2: Sex_M (M=1, F=0, M=1)
    REQUIRE(covariates(0, 2) == 1.0);
    REQUIRE(covariates(1, 2) == 0.0);
    REQUIRE(covariates(2, 2) == 1.0);
}
