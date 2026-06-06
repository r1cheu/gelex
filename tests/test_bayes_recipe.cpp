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

#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/half_normal_prior.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/recipe_options.h"
#include "gelex/data/genotype/genotype.h"
#include "gelex/exception.h"
#include "gelex/types/constrained_value.h"
#include "gelex/types/constrained_vector.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_effect_type.h"
#include "genotype_fixture.h"

using gelex::GelexException;
using gelex::GeneticMode;
using gelex::OpenUnitInterval;
using gelex::ScaleMultiplier;
using gelex::Simplex;
using gelex::bayes::BayesRecipe;
using gelex::bayes::BayesRecipeOptions;
using gelex::bayes::BayesRecipeScheme;
using gelex::bayes::to_bayes_recipe_scheme;

namespace
{

auto make_genotype(Eigen::MatrixXd data) -> gelex::genotype::Genotype
{
    auto mean = data.colwise().mean().transpose().eval();
    auto stddev = Eigen::VectorXd::Ones(data.cols());
    return gelex::test::GenotypeBuilder::build(
        std::move(data), std::move(mean), std::move(stddev));
}

auto make_model() -> gelex::BayesModel
{
    std::vector<gelex::bayes::GeneticDesign> genetics;
    genetics.emplace_back(
        GeneticMode::A,
        make_genotype(Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}}));
    genetics.emplace_back(
        GeneticMode::D,
        make_genotype(Eigen::MatrixXd{{1.0, 0.0}, {0.0, 1.0}, {1.0, 2.0}}));

    return gelex::BayesModel{
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        gelex::FixedDesign::build(3),
        {},
        std::move(genetics)};
}

}  // namespace

TEST_CASE("BayesRecipe scheme parser accepts short names", "[bayes_recipe]")
{
    REQUIRE(to_bayes_recipe_scheme("RR") == BayesRecipeScheme::RR);
    REQUIRE(to_bayes_recipe_scheme("A") == BayesRecipeScheme::A);
    REQUIRE(to_bayes_recipe_scheme("B") == BayesRecipeScheme::B);
    REQUIRE(to_bayes_recipe_scheme("C") == BayesRecipeScheme::C);
    REQUIRE(to_bayes_recipe_scheme("R") == BayesRecipeScheme::R);
    REQUIRE(to_bayes_recipe_scheme("CD") == BayesRecipeScheme::CD);

    REQUIRE_THROWS_AS(to_bayes_recipe_scheme("BayesRR"), GelexException);
    REQUIRE_THROWS_AS(to_bayes_recipe_scheme("BayesA"), GelexException);
    REQUIRE_THROWS_AS(to_bayes_recipe_scheme("BayesCD"), GelexException);
}

TEST_CASE("BayesR construction succeeds with defaults", "[bayes_recipe]")
{
    BayesRecipeOptions config;
    config.scheme = BayesRecipeScheme::R;
    config.modes = {GeneticMode::A};
    REQUIRE_NOTHROW(BayesRecipe(config));
}

TEST_CASE(
    "BayesRecipe rejects independent scheme-incompatible overrides",
    "[bayes_recipe]")
{
    SECTION("joint proportion")
    {
        BayesRecipeOptions config;
        config.scheme = BayesRecipeScheme::B;
        config.modes = {GeneticMode::A};
        config.joint_proportion = Simplex<double>{{0.5, 0.5}};
        REQUIRE_THROWS_AS(BayesRecipe(config), GelexException);
    }

    SECTION("dominance positive probability")
    {
        BayesRecipeOptions config;
        config.scheme = BayesRecipeScheme::B;
        config.modes = {GeneticMode::D};
        config.dominance_positive_probability
            = OpenUnitInterval<double>{0.6};
        REQUIRE_THROWS_AS(BayesRecipe(config), GelexException);
    }

    SECTION("proportion override")
    {
        BayesRecipeOptions config;
        config.modes = {GeneticMode::A};
        config.additive_proportion = Simplex<double>{{0.99, 0.01}};
        REQUIRE_THROWS_AS(BayesRecipe(config), GelexException);
    }

    SECTION("multiplier override")
    {
        BayesRecipeOptions config;
        config.scheme = BayesRecipeScheme::B;
        config.modes = {GeneticMode::A};
        config.additive_multiplier = ScaleMultiplier<double>{{0.0, 1.0}};
        REQUIRE_THROWS_AS(BayesRecipe(config), GelexException);
    }

    SECTION("unpaired BayesR proportion and multiplier")
    {
        BayesRecipeOptions config;
        config.scheme = BayesRecipeScheme::R;
        config.modes = {GeneticMode::A};
        config.additive_proportion = Simplex<double>{{0.99, 0.01}};
        REQUIRE_THROWS_AS(BayesRecipe(config), GelexException);
    }
}

TEST_CASE("BayesCD rejects joint_proportion of wrong size", "[bayes_recipe]")
{
    SECTION("2-element simplex rejected")
    {
        BayesRecipeOptions config;
        config.scheme = BayesRecipeScheme::CD;
        config.modes = {GeneticMode::A, GeneticMode::D};
        config.joint_proportion = Simplex<double>{{0.5, 0.5}};
        REQUIRE_THROWS_AS(BayesRecipe(config), GelexException);
    }
    SECTION("3-element simplex rejected")
    {
        BayesRecipeOptions config;
        config.scheme = BayesRecipeScheme::CD;
        config.modes = {GeneticMode::A, GeneticMode::D};
        config.joint_proportion = Simplex<double>{{0.8, 0.1, 0.1}};
        REQUIRE_THROWS_AS(BayesRecipe(config), GelexException);
    }
    SECTION("5-element simplex rejected")
    {
        BayesRecipeOptions config;
        config.scheme = BayesRecipeScheme::CD;
        config.modes = {GeneticMode::A, GeneticMode::D};
        config.joint_proportion
            = Simplex<double>{{0.8, 0.05, 0.05, 0.05, 0.05}};
        REQUIRE_THROWS_AS(BayesRecipe(config), GelexException);
    }
    SECTION("4-element simplex accepted")
    {
        BayesRecipeOptions config;
        config.scheme = BayesRecipeScheme::CD;
        config.modes = {GeneticMode::A, GeneticMode::D};
        config.joint_proportion = Simplex<double>{{0.9, 0.04, 0.03, 0.03}};
        REQUIRE_NOTHROW(BayesRecipe(config));
    }
}

TEST_CASE(
    "BayesRecipe make_prior creates scheme-specific genetic priors",
    "[bayes_recipe]")
{
    auto model = make_model();

    SECTION("RR creates a shared gaussian prior")
    {
        BayesRecipeOptions config;
        config.modes = {GeneticMode::A};
        auto prior = BayesRecipe(config).make_prior(model);
        auto genetics = prior.genetics();

        REQUIRE(genetics.size() == 1);
        const auto& single
            = std::get<gelex::bayes::SingleGeneticPrior>(genetics[0]);
        REQUIRE(
            std::holds_alternative<gelex::bayes::SingleSharedGaussianPrior>(
                single));
    }

    SECTION("A creates a per-marker gaussian prior")
    {
        BayesRecipeOptions config;
        config.scheme = BayesRecipeScheme::A;
        config.modes = {GeneticMode::A};
        auto prior = BayesRecipe(config).make_prior(model);
        auto genetics = prior.genetics();

        REQUIRE(genetics.size() == 1);
        const auto& single
            = std::get<gelex::bayes::SingleGeneticPrior>(genetics[0]);
        REQUIRE(
            std::holds_alternative<gelex::bayes::SinglePerMarkerGaussianPrior>(
                single));
    }

    SECTION("B creates a per-marker spike-slab gaussian prior")
    {
        BayesRecipeOptions config;
        config.scheme = BayesRecipeScheme::B;
        config.modes = {GeneticMode::A};
        auto prior = BayesRecipe(config).make_prior(model);
        auto genetics = prior.genetics();

        REQUIRE(genetics.size() == 1);
        const auto& single
            = std::get<gelex::bayes::SingleGeneticPrior>(genetics[0]);
        REQUIRE(
            std::holds_alternative<
                gelex::bayes::SinglePerMarkerSpikeSlabGaussianPrior>(single));
    }

    SECTION("C creates a shared spike-slab gaussian prior")
    {
        BayesRecipeOptions config;
        config.scheme = BayesRecipeScheme::C;
        config.modes = {GeneticMode::A};
        auto prior = BayesRecipe(config).make_prior(model);
        auto genetics = prior.genetics();

        REQUIRE(genetics.size() == 1);
        const auto& single
            = std::get<gelex::bayes::SingleGeneticPrior>(genetics[0]);
        REQUIRE(
            std::holds_alternative<
                gelex::bayes::SingleSharedSpikeSlabGaussianPrior>(single));
    }

    SECTION("R creates a scaled mixture gaussian prior")
    {
        BayesRecipeOptions config;
        config.scheme = BayesRecipeScheme::R;
        config.modes = {GeneticMode::A};
        auto prior = BayesRecipe(config).make_prior(model);
        auto genetics = prior.genetics();

        REQUIRE(genetics.size() == 1);
        const auto& single
            = std::get<gelex::bayes::SingleGeneticPrior>(genetics[0]);
        REQUIRE(
            std::holds_alternative<
                gelex::bayes::SingleScaledMixtureGaussianPrior>(single));
    }

    SECTION("CD creates a joint half normal mixture prior")
    {
        BayesRecipeOptions config;
        config.scheme = BayesRecipeScheme::CD;
        config.modes = {GeneticMode::A, GeneticMode::D};
        config.dominance_positive_probability = OpenUnitInterval<double>{0.7};
        auto prior = BayesRecipe(config).make_prior(model);
        auto genetics = prior.genetics();

        REQUIRE(genetics.size() == 1);
        const auto& joint
            = std::get<gelex::bayes::JointGeneticPrior>(genetics[0]);
        const auto& cd
            = std::get<gelex::bayes::JointHalfNormalMixturePrior>(joint);
        REQUIRE(cd.proportion().initial_value().isApprox(
            Eigen::VectorXd{
                {0.99,
                 0.003333333333333333,
                 0.003333333333333333,
                 0.003333333333333333}}));
        REQUIRE(cd.dominance_positive_probability().initial_value() == 0.7);
        REQUIRE(cd.dominance_positive_probability().prior().alpha() == 1.0);
        REQUIRE(cd.dominance_positive_probability().prior().beta() == 1.0);
    }
}

TEST_CASE(
    "BayesRecipe make_prior preserves independent mode order",
    "[bayes_recipe]")
{
    auto model = make_model();
    BayesRecipeOptions config;
    config.modes = {GeneticMode::A, GeneticMode::D};

    auto prior = BayesRecipe(config).make_prior(model);
    auto genetics = prior.genetics();

    REQUIRE(genetics.size() == 2);
    const auto& additive
        = std::get<gelex::bayes::SingleGeneticPrior>(genetics[0]);
    const auto& dominance
        = std::get<gelex::bayes::SingleGeneticPrior>(genetics[1]);
    REQUIRE(gelex::bayes::mode(additive) == GeneticMode::A);
    REQUIRE(gelex::bayes::mode(dominance) == GeneticMode::D);
}

TEST_CASE(
    "BayesRecipe accepts valid overrides for present modes",
    "[bayes_recipe]")
{
    SECTION("additive override with mode A present")
    {
        BayesRecipeOptions config;
        config.modes = {GeneticMode::A};
        config.additive_heritability = OpenUnitInterval<double>{0.3};
        REQUIRE_NOTHROW(BayesRecipe(config));
    }
    SECTION("dominance override with mode D present")
    {
        BayesRecipeOptions config;
        config.modes = {GeneticMode::D};
        config.dominance_heritability = OpenUnitInterval<double>{0.3};
        REQUIRE_NOTHROW(BayesRecipe(config));
    }
    SECTION("both overrides with modes {A, D}")
    {
        BayesRecipeOptions config;
        config.modes = {GeneticMode::A, GeneticMode::D};
        config.additive_heritability = OpenUnitInterval<double>{0.3};
        config.dominance_heritability = OpenUnitInterval<double>{0.2};
        REQUIRE_NOTHROW(BayesRecipe(config));
    }
}
