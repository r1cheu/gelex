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
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/recipe_options.h"
#include "gelex/data/genotype/genotype.h"
#include "gelex/exception.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/constrained_value.h"
#include "gelex/types/constrained_vector.h"
#include "gelex/types/genetic_effect_type.h"
#include "genotype_fixture.h"

using gelex::GelexException;
using gelex::GeneticMode;
using gelex::OpenUnitInterval;
using gelex::ScaleMultiplier;
using gelex::Simplex;
using gelex::bayes::BayesRecipe;
using gelex::bayes::BayesRecipeConfig;
using gelex::bayes::BayesRecipeScheme;
using gelex::bayes::EffectConfig;
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
    BayesRecipeConfig config;
    config.modes = {GeneticMode::A};
    REQUIRE_NOTHROW(BayesRecipe(BayesRecipeScheme::R, config));
}

TEST_CASE(
    "BayesRecipe rejects independent scheme-incompatible overrides",
    "[bayes_recipe]")
{
    SECTION("joint proportion")
    {
        BayesRecipeConfig config;
        config.modes = {GeneticMode::A};
        config.joint_proportion = Simplex<double>{{0.5, 0.5}};
        REQUIRE_THROWS_AS(
            BayesRecipe(BayesRecipeScheme::B, config), GelexException);
    }

    SECTION("proportion override")
    {
        BayesRecipeConfig config;
        config.modes = {GeneticMode::A};
        config.additive = EffectConfig{
            std::nullopt,
            std::optional<Simplex<double>>{Simplex<double>{{0.99, 0.01}}},
            std::nullopt,
            std::nullopt};
        REQUIRE_THROWS_AS(
            BayesRecipe(BayesRecipeScheme::RR, config), GelexException);
    }

    SECTION("multiplier override")
    {
        BayesRecipeConfig config;
        config.modes = {GeneticMode::A};
        config.additive = EffectConfig{
            std::nullopt,
            std::nullopt,
            std::optional<ScaleMultiplier<double>>{
                ScaleMultiplier<double>{{0.0, 1.0}}},
            std::nullopt};
        REQUIRE_THROWS_AS(
            BayesRecipe(BayesRecipeScheme::B, config), GelexException);
    }

    SECTION("unpaired BayesR proportion and multiplier")
    {
        BayesRecipeConfig config;
        config.modes = {GeneticMode::A};
        config.additive = EffectConfig{
            std::nullopt,
            std::optional<Simplex<double>>{Simplex<double>{{0.99, 0.01}}},
            std::nullopt,
            std::nullopt};
        REQUIRE_THROWS_AS(
            BayesRecipe(BayesRecipeScheme::R, config), GelexException);
    }
}

TEST_CASE("BayesCD rejects joint_proportion of wrong size", "[bayes_recipe]")
{
    SECTION("2-element simplex rejected")
    {
        BayesRecipeConfig config;
        config.modes = {GeneticMode::A, GeneticMode::D};
        config.joint_proportion = Simplex<double>{{0.5, 0.5}};
        REQUIRE_THROWS_AS(
            BayesRecipe(BayesRecipeScheme::CD, config), GelexException);
    }
    SECTION("3-element simplex rejected")
    {
        BayesRecipeConfig config;
        config.modes = {GeneticMode::A, GeneticMode::D};
        config.joint_proportion = Simplex<double>{{0.8, 0.1, 0.1}};
        REQUIRE_THROWS_AS(
            BayesRecipe(BayesRecipeScheme::CD, config), GelexException);
    }
    SECTION("5-element simplex rejected")
    {
        BayesRecipeConfig config;
        config.modes = {GeneticMode::A, GeneticMode::D};
        config.joint_proportion
            = Simplex<double>{{0.8, 0.05, 0.05, 0.05, 0.05}};
        REQUIRE_THROWS_AS(
            BayesRecipe(BayesRecipeScheme::CD, config), GelexException);
    }
    SECTION("4-element simplex accepted")
    {
        BayesRecipeConfig config;
        config.modes = {GeneticMode::A, GeneticMode::D};
        config.joint_proportion = Simplex<double>{{0.9, 0.04, 0.03, 0.03}};
        REQUIRE_NOTHROW(BayesRecipe(BayesRecipeScheme::CD, config));
    }
}

TEST_CASE(
    "BayesRecipe make_prior creates scheme-specific genetic priors",
    "[bayes_recipe]")
{
    auto model = make_model();

    SECTION("RR creates a shared gaussian prior")
    {
        BayesRecipeConfig config;
        config.modes = {GeneticMode::A};
        auto prior = BayesRecipe(BayesRecipeScheme::RR, config).make_prior(
            model);
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
        BayesRecipeConfig config;
        config.modes = {GeneticMode::A};
        auto prior = BayesRecipe(BayesRecipeScheme::A, config).make_prior(
            model);
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
        BayesRecipeConfig config;
        config.modes = {GeneticMode::A};
        auto prior = BayesRecipe(BayesRecipeScheme::B, config).make_prior(
            model);
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
        BayesRecipeConfig config;
        config.modes = {GeneticMode::A};
        auto prior = BayesRecipe(BayesRecipeScheme::C, config).make_prior(
            model);
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
        BayesRecipeConfig config;
        config.modes = {GeneticMode::A};
        auto prior = BayesRecipe(BayesRecipeScheme::R, config).make_prior(
            model);
        auto genetics = prior.genetics();

        REQUIRE(genetics.size() == 1);
        const auto& single
            = std::get<gelex::bayes::SingleGeneticPrior>(genetics[0]);
        REQUIRE(
            std::holds_alternative<
                gelex::bayes::SingleScaledMixtureGaussianPrior>(single));
    }

    SECTION("CD creates a joint gaussian mixture prior")
    {
        BayesRecipeConfig config;
        config.modes = {GeneticMode::A, GeneticMode::D};
        auto prior = BayesRecipe(BayesRecipeScheme::CD, config).make_prior(
            model);
        auto genetics = prior.genetics();

        REQUIRE(genetics.size() == 1);
        const auto& joint
            = std::get<gelex::bayes::JointGeneticPrior>(genetics[0]);
        REQUIRE(std::holds_alternative<gelex::bayes::JointGaussianMixturePrior>(
            joint));
    }
}

TEST_CASE(
    "BayesRecipe make_prior preserves independent mode order",
    "[bayes_recipe]")
{
    auto model = make_model();
    BayesRecipeConfig config;
    config.modes = {GeneticMode::A, GeneticMode::D};

    auto prior = BayesRecipe(BayesRecipeScheme::RR, config).make_prior(model);
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
        BayesRecipeConfig config;
        config.modes = {GeneticMode::A};
        config.additive = EffectConfig{
            std::optional<OpenUnitInterval<double>>{std::in_place, 0.3},
            std::nullopt,
            std::nullopt,
            std::nullopt};
        REQUIRE_NOTHROW(BayesRecipe(BayesRecipeScheme::RR, config));
    }
    SECTION("dominance override with mode D present")
    {
        BayesRecipeConfig config;
        config.modes = {GeneticMode::D};
        config.dominance = EffectConfig{
            std::optional<OpenUnitInterval<double>>{std::in_place, 0.3},
            std::nullopt,
            std::nullopt,
            std::nullopt};
        REQUIRE_NOTHROW(BayesRecipe(BayesRecipeScheme::RR, config));
    }
    SECTION("both overrides with modes {A, D}")
    {
        BayesRecipeConfig config;
        config.modes = {GeneticMode::A, GeneticMode::D};
        config.additive = EffectConfig{
            std::optional<OpenUnitInterval<double>>{std::in_place, 0.3},
            std::nullopt,
            std::nullopt,
            std::nullopt};
        config.dominance = EffectConfig{
            std::optional<OpenUnitInterval<double>>{std::in_place, 0.2},
            std::nullopt,
            std::nullopt,
            std::nullopt};
        REQUIRE_NOTHROW(BayesRecipe(BayesRecipeScheme::RR, config));
    }
}
