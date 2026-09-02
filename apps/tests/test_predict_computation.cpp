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
#include <vector>

#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

#include "cli/predict/compute.h"

using gelex::GeneticMode;

// ---------------------------------------------------------------------------
// compute_gebv tests
// ---------------------------------------------------------------------------

TEST_CASE("compute_gebv: additive only", "[predict][compute]")
{
    // add = [[1,2],[3,4]], effects.add = [0.5, -0.5]
    // add_predictions = [1*0.5 + 2*(-0.5), 3*0.5 + 4*(-0.5)] = [-0.5, -0.5]
    gelex::ModeMap<Eigen::MatrixXd> geno;
    geno[GeneticMode::A] = Eigen::MatrixXd{{1.0, 2.0}, {3.0, 4.0}};

    gelex::ModeMap<Eigen::VectorXd> effects;
    effects[GeneticMode::A] = Eigen::VectorXd{{0.5, -0.5}};

    const auto result = cli::compute_gebv(geno, effects);

    Eigen::VectorXd expected_add{{-0.5, -0.5}};
    REQUIRE(result.at(GeneticMode::A).isApprox(expected_add));
    REQUIRE_FALSE(result.contains(GeneticMode::D));
}

TEST_CASE("compute_gebv: dominance only", "[predict][compute]")
{
    // dom = [[1,0],[0,1]], effects.dom = [1.0, 2.0]
    // dom_predictions = [1.0, 2.0]; no additive component
    gelex::ModeMap<Eigen::MatrixXd> geno;
    geno[GeneticMode::D] = Eigen::MatrixXd{{1.0, 0.0}, {0.0, 1.0}};

    gelex::ModeMap<Eigen::VectorXd> effects;
    effects[GeneticMode::D] = Eigen::VectorXd{{1.0, 2.0}};

    const auto result = cli::compute_gebv(geno, effects);

    Eigen::VectorXd expected_dom{{1.0, 2.0}};
    REQUIRE_FALSE(result.contains(GeneticMode::A));
    REQUIRE(result.at(GeneticMode::D).isApprox(expected_dom));
}

TEST_CASE("compute_gebv: additive + dominance", "[predict][compute]")
{
    // add = [[1,0],[0,1]], effects.add = [1.0, 2.0]
    // dom = [[0.5,0.5],[0.5,0.5]], effects.dom = [0.1, 0.2]
    // add_pred = [1*1.0+0*2.0, 0*1.0+1*2.0] = [1.0, 2.0]
    // dom_pred = [0.5*0.1+0.5*0.2, 0.5*0.1+0.5*0.2] = [0.15, 0.15]
    // total = [1.15, 2.15]
    gelex::ModeMap<Eigen::MatrixXd> geno;
    geno[GeneticMode::A] = Eigen::MatrixXd{{1.0, 0.0}, {0.0, 1.0}};
    geno[GeneticMode::D] = Eigen::MatrixXd{{0.5, 0.5}, {0.5, 0.5}};

    gelex::ModeMap<Eigen::VectorXd> effects;
    effects[GeneticMode::A] = Eigen::VectorXd{{1.0, 2.0}};
    effects[GeneticMode::D] = Eigen::VectorXd{{0.1, 0.2}};

    const auto result = cli::compute_gebv(geno, effects);

    Eigen::VectorXd expected_add{{1.0, 2.0}};
    Eigen::VectorXd expected_dom{{0.15, 0.15}};
    REQUIRE(result.at(GeneticMode::A).isApprox(expected_add));
    REQUIRE(result.at(GeneticMode::D).isApprox(expected_dom));
}

// ---------------------------------------------------------------------------
// compute_covariate_effects tests
// ---------------------------------------------------------------------------

TEST_CASE("compute_covariate_effects: intercept only", "[predict][compute]")
{
    // covariates = [[1],[1],[1]], coefficients = {["Intercept"], [2.5]}
    // intercept is excluded: per_covariate has 0 cols, covar_names empty
    Eigen::MatrixXd covariates{{1.0}, {1.0}, {1.0}};

    std::vector<std::string> term_names{"Intercept"};
    Eigen::VectorXd coefficients{{2.5}};

    const auto result
        = cli::compute_covariate_effects(covariates, term_names, coefficients);

    REQUIRE(result.per_covariate.cols() == 0);
    REQUIRE(result.covar_names.empty());
}

TEST_CASE(
    "compute_covariate_effects: intercept + covariates",
    "[predict][compute]")
{
    // covariates = [[1, 25, 1], [1, 30, 0]]
    // coefficients = {["Intercept","Age","Sex_M"], [1.0, 0.2, -0.3]}
    // per_covariate: Age=[5.0, 6.0], Sex_M=[-0.3, 0.0]
    Eigen::MatrixXd covariates{{1.0, 25.0, 1.0}, {1.0, 30.0, 0.0}};

    std::vector<std::string> term_names{"Intercept", "Age", "Sex_M"};
    Eigen::VectorXd coefficients{{1.0, 0.2, -0.3}};

    const auto result
        = cli::compute_covariate_effects(covariates, term_names, coefficients);

    REQUIRE(result.covar_names == std::vector<std::string>{"Age", "Sex_M"});
    REQUIRE(result.per_covariate.cols() == 2);

    Eigen::VectorXd expected_age{{5.0, 6.0}};
    Eigen::VectorXd expected_sex{{-0.3, 0.0}};
    REQUIRE(result.per_covariate.col(0).isApprox(expected_age));
    REQUIRE(result.per_covariate.col(1).isApprox(expected_sex));
}

TEST_CASE(
    "compute_covariate_effects: rejects terms not led by Intercept",
    "[predict][compute]")
{
    Eigen::MatrixXd covariates{{25.0, 1.0}, {30.0, 1.0}};

    std::vector<std::string> term_names{"Age", "Intercept"};
    Eigen::VectorXd coefficients{{0.2, 1.0}};

    REQUIRE_THROWS_AS(
        cli::compute_covariate_effects(covariates, term_names, coefficients),
        gelex::GelexException);
}
