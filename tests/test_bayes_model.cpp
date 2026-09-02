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
#include <utility>
#include <vector>

#include "gelex/bayes/model.h"
#include "gelex/data/genotype_method.h"
#include "gelex/exception.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_mode.h"

#include "compact_genotype_fixture.h"

using gelex::BayesModel;
using gelex::FixedDesign;
using gelex::GelexException;
using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::GenotypeMethod;

namespace
{

const Eigen::MatrixXd GENOTYPES{{0.0, 0.0}, {1.0, 1.0}, {2.0, 1.0}};

auto make_phenotype() -> Eigen::VectorXd
{
    return Eigen::VectorXd{{1.0, 2.0, 4.0}};
}

}  // namespace

TEST_CASE("BayesModel rejects design row mismatches", "[bayes_model]")
{
    SECTION("fixed")
    {
        REQUIRE_THROWS_AS(
            BayesModel(
                make_phenotype(),
                FixedDesign::make(2),
                {},
                gelex::test::make_genetic_design(GENOTYPES)),
            GelexException);
    }

    SECTION("random")
    {
        std::vector<gelex::bayes::RandomDesign> random;
        random.emplace_back(
            "batch",
            std::vector<std::string>{"a", "b"},
            Eigen::MatrixXd{{1.0, 0.0}, {0.0, 1.0}});

        REQUIRE_THROWS_AS(
            BayesModel(
                make_phenotype(),
                FixedDesign::make(3),
                std::move(random),
                gelex::test::make_genetic_design(GENOTYPES)),
            GelexException);
    }

    SECTION("genetic")
    {
        REQUIRE_THROWS_AS(
            BayesModel(
                make_phenotype(),
                FixedDesign::make(3),
                {},
                gelex::test::make_genetic_design(
                    Eigen::MatrixXd{{0.0}, {1.0}})),
            GelexException);
    }
}

TEST_CASE("BayesModel accepts a design without projections", "[bayes_model]")
{
    auto model = BayesModel{
        make_phenotype(),
        FixedDesign::make(3),
        {},
        gelex::test::make_genetic_design(GENOTYPES, GeneticModeSet{})};

    REQUIRE(model.genetic().modes().size() == 0);
}

TEST_CASE("GeneticDesign exposes compact column metadata", "[bayes_model]")
{
    auto model = gelex::test::make_compact_model(
        GENOTYPES,
        make_phenotype(),
        GeneticModeSet{GeneticMode::A},
        GenotypeMethod::Center);
    const auto& design = model.genetic();
    const Eigen::MatrixXd expected{
        {-1.0, -2.0 / 3.0}, {0.0, 1.0 / 3.0}, {1.0, 1.0 / 3.0}};

    REQUIRE(design.xtx_diag().isApprox(
        expected.colwise().squaredNorm().transpose()));
    REQUIRE(
        design.col_var().isApprox(Eigen::RowVectorXd{{2.0 / 3.0, 2.0 / 9.0}}));
    REQUIRE(design.valid_indices().size() == 2);
}

TEST_CASE("Moved BayesModel keeps compact design valid", "[bayes_model]")
{
    auto source = gelex::test::make_compact_model(
        GENOTYPES,
        make_phenotype(),
        GeneticMode::A | GeneticMode::D,
        GenotypeMethod::OrthCenter);
    auto model = std::move(source);
    const auto& design = model.genetic();
    const Eigen::VectorXd values{{1.0, -0.5, 0.25}};
    Eigen::VectorXd additive_output = Eigen::VectorXd::Zero(3);
    Eigen::VectorXd dominance_output = Eigen::VectorXd::Zero(3);

    design.axpy(GeneticMode::A, 0, 1.0, additive_output);
    design.axpy(GeneticMode::D, 0, 1.0, dominance_output);

    REQUIRE(
        design.dot(GeneticMode::A, 0, values) == additive_output.dot(values));
    REQUIRE(
        design.dot(GeneticMode::D, 0, values) == dominance_output.dot(values));
    REQUIRE(design.modes().size() == 2);
}
