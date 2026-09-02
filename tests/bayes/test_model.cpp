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
#include <catch2/matchers/catch_matchers_string.hpp>
#include <limits>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

#include "gelex/bayes/model.h"
#include "gelex/data/fixed_design.h"
#include "gelex/data/genotype_method.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

#include "compact_genotype_fixture.h"
#include "random_design_fixture.h"

using Catch::Matchers::ContainsSubstring;
using gelex::BayesModel;
using gelex::FixedDesign;
using gelex::GelexException;
using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::GenotypeMethod;

namespace
{

const Eigen::MatrixXd genotypes{{0.0, 0.0}, {1.0, 1.0}, {2.0, 1.0}};

auto make_phenotype() -> Eigen::VectorXd
{
    return Eigen::VectorXd{{1.0, 2.0, 4.0}};
}

auto make_model(Eigen::VectorXd phenotype) -> BayesModel
{
    return BayesModel{
        std::move(phenotype),
        FixedDesign::make(genotypes.rows()),
        {},
        gelex::test::make_genetic_design(genotypes)};
}

}  // namespace

TEST_CASE(
    "BayesModel requires a finite positive phenotype variance",
    "[bayes_model]")
{
    SECTION("empty phenotype")
    {
        REQUIRE_THROWS_WITH(
            make_model(Eigen::VectorXd{}),
            ContainsSubstring("phenotype must not be empty"));
    }

    SECTION("constant phenotype")
    {
        REQUIRE_THROWS_WITH(
            make_model(Eigen::VectorXd{{1.0, 1.0, 1.0}}),
            ContainsSubstring(
                "phenotype variance must be finite and positive"));
    }

    SECTION("non-finite phenotype")
    {
        REQUIRE_THROWS_WITH(
            make_model(
                Eigen::VectorXd{
                    {1.0, std::numeric_limits<double>::quiet_NaN(), 3.0}}),
            ContainsSubstring(
                "phenotype variance must be finite and positive"));
    }
}

TEST_CASE("BayesModel rejects design row mismatches", "[bayes_model]")
{
    SECTION("fixed")
    {
        REQUIRE_THROWS_AS(
            BayesModel(
                make_phenotype(),
                FixedDesign::make(2),
                {},
                gelex::test::make_genetic_design(genotypes)),
            GelexException);
    }

    SECTION("random")
    {
        std::vector<gelex::bayes::RandomDesign> random;
        random.push_back(
            gelex::test::make_random_design(
                "batch",
                std::vector<std::string>{"a", "b"},
                Eigen::MatrixXd{{1.0, 0.0}, {0.0, 1.0}}));

        REQUIRE_THROWS_AS(
            BayesModel(
                make_phenotype(),
                FixedDesign::make(3),
                std::move(random),
                gelex::test::make_genetic_design(genotypes)),
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
        gelex::test::make_genetic_design_without_modes(genotypes)};

    REQUIRE(std::ranges::empty(model.genetic().each_mode()));
}

TEST_CASE("GeneticDesign exposes compact column metadata", "[bayes_model]")
{
    auto model = gelex::test::make_compact_model(
        genotypes,
        make_phenotype(),
        GeneticModeSet{GeneticMode::A},
        GenotypeMethod::Center);
    const auto& design = model.genetic();
    const auto& projection = design.projection(GeneticMode::A);
    const Eigen::MatrixXd expected{
        {-1.0, -2.0 / 3.0}, {0.0, 1.0 / 3.0}, {1.0, 1.0 / 3.0}};

    REQUIRE(projection.xtx_diag().isApprox(
        expected.colwise().squaredNorm().transpose()));
    REQUIRE(projection.col_var().isApprox(
        Eigen::RowVectorXd{{2.0 / 3.0, 2.0 / 9.0}}));
    REQUIRE(projection.valid_indices().size() == 2);
}

TEST_CASE("Moved BayesModel keeps compact design valid", "[bayes_model]")
{
    auto source = gelex::test::make_compact_model(
        genotypes,
        make_phenotype(),
        GeneticMode::A | GeneticMode::D,
        GenotypeMethod::OrthCenter);
    auto model = std::move(source);
    const auto& design = model.genetic();
    const Eigen::VectorXd values{{1.0, -0.5, 0.25}};
    Eigen::VectorXd additive_output = Eigen::VectorXd::Zero(3);
    Eigen::VectorXd dominance_output = Eigen::VectorXd::Zero(3);
    const auto& additive_projection = design.projection(GeneticMode::A);
    const auto& dominance_projection = design.projection(GeneticMode::D);

    additive_projection.axpy(0, 1.0, additive_output);
    dominance_projection.axpy(0, 1.0, dominance_output);

    REQUIRE(additive_projection.dot(0, values) == additive_output.dot(values));
    REQUIRE(
        dominance_projection.dot(0, values) == dominance_output.dot(values));
    REQUIRE(design.contains(GeneticMode::A));
    REQUIRE(design.contains(GeneticMode::D));
}
