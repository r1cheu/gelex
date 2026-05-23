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

#include <catch2/catch_test_macros.hpp>

#include "gelex/exception.h"
#include "gelex/model/bayes/recipe.h"
#include "gelex/model/bayes/recipe_options.h"
#include "gelex/types/constrained_value.h"
#include "gelex/types/constrained_vector.h"
#include "gelex/types/genetic_effect_type.h"

using gelex::GelexException;
using gelex::GeneticMode;
using gelex::OpenUnitInterval;
using gelex::Simplex;
using gelex::bayes::BayesRecipe;
using gelex::bayes::BayesRecipeConfig;
using gelex::bayes::BayesRecipePreset;
using gelex::bayes::EffectConfig;
using gelex::bayes::to_bayes_recipe_preset;

TEST_CASE("BayesRecipe preset parser accepts short names", "[bayes_recipe]")
{
    REQUIRE(to_bayes_recipe_preset("RR") == BayesRecipePreset::RR);
    REQUIRE(to_bayes_recipe_preset("A") == BayesRecipePreset::A);
    REQUIRE(to_bayes_recipe_preset("B") == BayesRecipePreset::B);
    REQUIRE(to_bayes_recipe_preset("C") == BayesRecipePreset::C);
    REQUIRE(to_bayes_recipe_preset("R") == BayesRecipePreset::R);
    REQUIRE(to_bayes_recipe_preset("CD") == BayesRecipePreset::CD);

    REQUIRE_THROWS_AS(to_bayes_recipe_preset("BayesRR"), GelexException);
    REQUIRE_THROWS_AS(to_bayes_recipe_preset("BayesA"), GelexException);
    REQUIRE_THROWS_AS(to_bayes_recipe_preset("BayesCD"), GelexException);
}

TEST_CASE("BayesR construction succeeds with defaults", "[bayes_recipe]")
{
    BayesRecipeConfig config;
    config.modes = {GeneticMode::A};
    REQUIRE_NOTHROW(BayesRecipe(BayesRecipePreset::R, config));
}

TEST_CASE("BayesCD rejects joint_proportion of wrong size", "[bayes_recipe]")
{
    SECTION("2-element simplex rejected")
    {
        BayesRecipeConfig config;
        config.modes = {GeneticMode::A, GeneticMode::D};
        config.joint_proportion = Simplex<double>{{0.5, 0.5}};
        REQUIRE_THROWS_AS(
            BayesRecipe(BayesRecipePreset::CD, config), GelexException);
    }
    SECTION("3-element simplex rejected")
    {
        BayesRecipeConfig config;
        config.modes = {GeneticMode::A, GeneticMode::D};
        config.joint_proportion = Simplex<double>{{0.8, 0.1, 0.1}};
        REQUIRE_THROWS_AS(
            BayesRecipe(BayesRecipePreset::CD, config), GelexException);
    }
    SECTION("5-element simplex rejected")
    {
        BayesRecipeConfig config;
        config.modes = {GeneticMode::A, GeneticMode::D};
        config.joint_proportion
            = Simplex<double>{{0.8, 0.05, 0.05, 0.05, 0.05}};
        REQUIRE_THROWS_AS(
            BayesRecipe(BayesRecipePreset::CD, config), GelexException);
    }
    SECTION("4-element simplex accepted")
    {
        BayesRecipeConfig config;
        config.modes = {GeneticMode::A, GeneticMode::D};
        config.joint_proportion
            = Simplex<double>{{0.9, 0.04, 0.03, 0.03}};
        REQUIRE_NOTHROW(BayesRecipe(BayesRecipePreset::CD, config));
    }
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
        REQUIRE_NOTHROW(BayesRecipe(BayesRecipePreset::RR, config));
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
        REQUIRE_NOTHROW(BayesRecipe(BayesRecipePreset::RR, config));
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
        REQUIRE_NOTHROW(BayesRecipe(BayesRecipePreset::RR, config));
    }
}
