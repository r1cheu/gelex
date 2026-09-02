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

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <concepts>
#include <type_traits>
#include <utility>

#include "gelex/bayes/builtin_method.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

#include "cli/mcmc/config.h"
#include "cli/mcmc/recipe.h"

namespace
{

constexpr auto mode_a = gelex::GeneticModeSet{gelex::GeneticMode::A};
constexpr auto mode_d = gelex::GeneticModeSet{gelex::GeneticMode::D};
constexpr auto mode_ad = gelex::GeneticMode::A | gelex::GeneticMode::D;

template <
    gelex::GeneticModeSet Modes,
    gelex::BayesMethod Method,
    typename Inspect>
auto inspect_mcmc_recipe(cli::McmcConfig config, Inspect inspect) -> void
{
    config.mode = Modes;
    config.method = Method;
    auto matched = false;
    cli::dispatch_mcmc_recipe(
        config,
        [&](const auto& recipe)
        {
            using Recipe = std::remove_cvref_t<decltype(recipe)>;
            using Expected = gelex::BuiltinBayesRecipe<Modes, Method>;
            if constexpr (std::same_as<Recipe, Expected>)
            {
                matched = true;
                inspect(recipe);
            }
        });
    REQUIRE(matched);
}

}  // namespace

TEST_CASE(
    "MCMC recipe adapter maps variance and spike-slab options",
    "[cli][mcmc][recipe]")
{
    cli::McmcConfig config;
    config.genetic_variance_shares.get<gelex::GeneticMode::A>() = 0.3;
    config.genetic_variance_shares.get<gelex::GeneticMode::D>() = 0.1;
    config.mixture_probabilities.get<gelex::GeneticMode::A>() = {0.2};
    config.mixture_probabilities.get<gelex::GeneticMode::D>() = {0.4};

    inspect_mcmc_recipe<mode_ad, gelex::BayesMethod::B>(
        std::move(config),
        [](const auto& recipe)
        {
            REQUIRE(recipe.variance().share(gelex::GeneticMode::A) == 0.3);
            REQUIRE(recipe.variance().share(gelex::GeneticMode::D) == 0.1);
            REQUIRE(recipe.variance().residual() == 0.6);
            REQUIRE(
                recipe.genetic_spec()
                    .template get<gelex::GeneticMode::A>()
                    .probability()
                == 0.2);
            REQUIRE(
                recipe.genetic_spec()
                    .template get<gelex::GeneticMode::D>()
                    .probability()
                == 0.4);
        });
}

TEST_CASE(
    "MCMC recipe adapter maps scaled-mixture and joint options",
    "[cli][mcmc][recipe]")
{
    cli::McmcConfig scaled_config;
    scaled_config.mixture_probabilities.get<gelex::GeneticMode::A>()
        = {0.5, 0.2, 0.15, 0.1, 0.05};
    scaled_config.mixture_scales.get<gelex::GeneticMode::A>()
        = {0.0, 0.001, 0.01, 0.1, 1.0};
    inspect_mcmc_recipe<mode_a, gelex::BayesMethod::R>(
        std::move(scaled_config),
        [](const auto& recipe)
        {
            REQUIRE(
                recipe.genetic_spec()
                    .template get<gelex::GeneticMode::A>()
                    .probabilities()
                == std::array{0.5, 0.2, 0.15, 0.1, 0.05});
            REQUIRE(
                recipe.genetic_spec()
                    .template get<gelex::GeneticMode::A>()
                    .scales()
                == std::array{0.0, 0.001, 0.01, 0.1, 1.0});
        });

    cli::McmcConfig joint_config;
    joint_config.dominance_positive_probability = 0.7;
    joint_config.mixture_probabilities.joint() = {0.7, 0.1, 0.1, 0.1};
    inspect_mcmc_recipe<mode_ad, gelex::BayesMethod::CD>(
        std::move(joint_config),
        [](const auto& recipe)
        {
            REQUIRE(
                recipe.genetic_spec()
                    .template get<gelex::GeneticMode::D>()
                    .positive_probability()
                == 0.7);
            REQUIRE(
                recipe.genetic_spec().joint().probabilities()
                == std::array{0.7, 0.1, 0.1, 0.1});
        });
}

TEST_CASE(
    "MCMC recipe dispatch selects every supported method and mode",
    "[cli][mcmc][recipe]")
{
    constexpr std::array modes{mode_a, mode_d, mode_ad};
    for (const auto mode : modes)
    {
        for (const auto [method, name] : gelex::bayes_method_names)
        {
            static_cast<void>(name);
            cli::McmcConfig config;
            config.mode = mode;
            config.method = method;
            if (method == gelex::BayesMethod::CD && mode != mode_ad)
            {
                REQUIRE_THROWS_AS(
                    cli::dispatch_mcmc_recipe(config, [](const auto&) {}),
                    gelex::GelexException);
            }
            else
            {
                REQUIRE_NOTHROW(
                    cli::dispatch_mcmc_recipe(config, [](const auto&) {}));
            }
        }
    }
}

TEST_CASE(
    "MCMC config rejects values for absent genetic modes",
    "[cli][mcmc][recipe]")
{
    cli::McmcConfig config;

    SECTION("variance share")
    {
        config.mode = mode_d;
        config.genetic_variance_shares.get<gelex::GeneticMode::A>() = 0.2;
    }
    SECTION("mixture probabilities")
    {
        config.method = gelex::BayesMethod::B;
        config.mixture_probabilities.get<gelex::GeneticMode::D>() = {0.2};
    }
    SECTION("mixture scales")
    {
        config.method = gelex::BayesMethod::R;
        config.mixture_scales.get<gelex::GeneticMode::D>()
            = {0.0, 0.001, 0.01, 0.1, 1.0};
    }
    SECTION("dominance positive probability")
    {
        config.dominance_positive_probability = 0.7;
    }

    REQUIRE_THROWS_AS(
        cli::dispatch_mcmc_recipe(config, [](const auto&) {}),
        gelex::GelexException);
}

TEST_CASE(
    "MCMC config rejects genetic options unsupported by the method",
    "[cli][mcmc][recipe]")
{
    cli::McmcConfig config;
    config.mode = mode_ad;

    SECTION("mode probabilities")
    {
        config.mixture_probabilities.get<gelex::GeneticMode::A>() = {0.2};
    }
    SECTION("mode scales")
    {
        config.method = gelex::BayesMethod::B;
        config.mixture_scales.get<gelex::GeneticMode::A>()
            = {0.0, 0.001, 0.01, 0.1, 1.0};
    }
    SECTION("joint probabilities")
    {
        config.method = gelex::BayesMethod::R;
        config.mixture_probabilities.joint() = {0.7, 0.1, 0.1, 0.1};
    }
    SECTION("dominance positive probability")
    {
        config.method = gelex::BayesMethod::R;
        config.dominance_positive_probability = 0.7;
    }

    REQUIRE_THROWS_AS(
        cli::dispatch_mcmc_recipe(config, [](const auto&) {}),
        gelex::GelexException);
}

TEST_CASE(
    "MCMC recipe adapter validates family-specific vector sizes",
    "[cli][mcmc][recipe]")
{
    cli::McmcConfig config;

    SECTION("spike-slab probability")
    {
        config.method = gelex::BayesMethod::B;
        config.mixture_probabilities.get<gelex::GeneticMode::A>() = {0.1, 0.9};
    }
    SECTION("scaled-mixture probability")
    {
        config.method = gelex::BayesMethod::R;
        config.mixture_probabilities.get<gelex::GeneticMode::A>() = {0.2, 0.2};
    }
    SECTION("joint probability")
    {
        config.mode = mode_ad;
        config.method = gelex::BayesMethod::CD;
        config.mixture_probabilities.joint() = {0.5, 0.5};
    }

    REQUIRE_THROWS_AS(
        cli::dispatch_mcmc_recipe(config, [](const auto&) {}),
        gelex::GelexException);
}

TEST_CASE(
    "MCMC recipe adapter requires random inputs and variance together",
    "[cli][mcmc][recipe]")
{
    cli::McmcConfig config;
    config.random_pve = 0.1;
    REQUIRE_THROWS_AS(
        cli::dispatch_mcmc_recipe(config, [](const auto&) {}),
        gelex::GelexException);

    config.random.qrand_paths = {"random.tsv"};
    inspect_mcmc_recipe<mode_a, gelex::BayesMethod::RR>(
        config,
        [](const auto& recipe) { REQUIRE(recipe.variance().random() == 0.1); });

    config.random_pve.reset();
    REQUIRE_THROWS_AS(
        cli::dispatch_mcmc_recipe(config, [](const auto&) {}),
        gelex::GelexException);
}
