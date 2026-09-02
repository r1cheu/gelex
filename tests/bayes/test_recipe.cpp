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

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <concepts>
#include <string>
#include <type_traits>

#include "gelex/bayes/builtin_recipe.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/semantic_method.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/variance_budget.h"
#include "gelex/exception.h"
#include "gelex/types/genetic_mode.h"
#include "gelex/types/mode_values.h"

using Catch::Matchers::ContainsSubstring;
using gelex::BayesMethod;
using gelex::BayesRecipe;
using gelex::BuiltinBayesRecipe;
using gelex::Gaussian;
using gelex::GaussianMethod;
using gelex::GelexException;
using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::JointSpikeSlab;
using gelex::JointSpikeSlabMethod;
using gelex::ModeValues;
using gelex::ScaledMixture;
using gelex::ScaledMixtureMethod;
using gelex::SpikeSlab;
using gelex::SpikeSlabMethod;
using gelex::UpdatePolicy;
using gelex::VarianceBudget;
using gelex::VarianceLayout;

namespace
{

constexpr auto mode_a = GeneticModeSet{GeneticMode::A};
constexpr auto mode_ad = GeneticMode::A | GeneticMode::D;

using SpikeSlabA = ModeValues<mode_a, SpikeSlab>;
using SpikeSlabAD = ModeValues<mode_ad, SpikeSlab, SpikeSlab>;
using ScaledMixtureAD = ModeValues<mode_ad, ScaledMixture, ScaledMixture>;
using JointSpikeSlabAD = gelex::JointModeValues<
    ModeValues<
        mode_ad,
        Gaussian,
        gelex::HalfNormal<gelex::HalfNormalAsymmetry::Count>>,
    JointSpikeSlab>;
using PooledGaussianMethod = GaussianMethod<VarianceLayout::Pooled>;
using UnpooledGaussianMethod = GaussianMethod<VarianceLayout::Unpooled>;
using PooledSpikeSlabMethod = SpikeSlabMethod<VarianceLayout::Pooled>;
using UnpooledSpikeSlabMethod = SpikeSlabMethod<VarianceLayout::Unpooled>;
using FixedPooledSpikeSlabMethod
    = SpikeSlabMethod<VarianceLayout::Pooled, UpdatePolicy::Fixed>;
using DefaultScaledMixtureMethod = ScaledMixtureMethod<>;
using DefaultJointSpikeSlabMethod = JointSpikeSlabMethod<>;
using MagnitudeJointSpikeSlabMethod = JointSpikeSlabMethod<
    UpdatePolicy::Sampled,
    gelex::HalfNormalAsymmetry::Magnitude>;

template <GeneticModeSet Modes, typename Method>
concept RecipeExists = requires { typename BayesRecipe<Modes, Method>; };

template <typename Recipe>
concept HasMethodToken = requires { Recipe::method; };

static_assert(BayesRecipe<mode_a, DefaultScaledMixtureMethod>::modes == mode_a);
static_assert(!HasMethodToken<BayesRecipe<mode_a, DefaultScaledMixtureMethod>>);
static_assert(std::same_as<
              BayesRecipe<mode_a, PooledGaussianMethod>::method_type,
              PooledGaussianMethod>);
static_assert(std::same_as<
              BuiltinBayesRecipe<mode_a, BayesMethod::RR>,
              BayesRecipe<mode_a, PooledGaussianMethod>>);
static_assert(std::same_as<
              BuiltinBayesRecipe<mode_a, BayesMethod::A>,
              BayesRecipe<mode_a, UnpooledGaussianMethod>>);
static_assert(std::same_as<
              BuiltinBayesRecipe<mode_a, BayesMethod::B>,
              BayesRecipe<mode_a, UnpooledSpikeSlabMethod>>);
static_assert(std::same_as<
              BuiltinBayesRecipe<mode_a, BayesMethod::C>,
              BayesRecipe<mode_a, PooledSpikeSlabMethod>>);
static_assert(std::same_as<
              BuiltinBayesRecipe<mode_a, BayesMethod::R>,
              BayesRecipe<mode_a, DefaultScaledMixtureMethod>>);
static_assert(std::same_as<
              BuiltinBayesRecipe<mode_ad, BayesMethod::CD>,
              BayesRecipe<mode_ad, DefaultJointSpikeSlabMethod>>);

// B and C take identical genetic specs but stay distinct types: the shared and
// per-marker variance difference is intrinsic to the method.
static_assert(!std::same_as<
              BayesRecipe<mode_ad, UnpooledSpikeSlabMethod>,
              BayesRecipe<mode_ad, PooledSpikeSlabMethod>>);

// RR and A have empty genetic specs: their prior shape is fixed by the data
// alone. Gaussian says so at the type level and costs no storage.
static_assert(
    std::same_as<
        typename BayesRecipe<mode_ad, PooledGaussianMethod>::genetic_spec_type,
        Gaussian>);
static_assert(
    std::same_as<
        std::remove_cvref_t<
            decltype(BayesRecipe<mode_ad, UnpooledGaussianMethod>::defaults()
                         .genetic_spec())>,
        Gaussian>);
static_assert(
    sizeof(BayesRecipe<mode_ad, PooledGaussianMethod>)
    == sizeof(VarianceBudget));
static_assert(std::constructible_from<
              BayesRecipe<mode_ad, PooledGaussianMethod>,
              VarianceBudget>);
static_assert(!std::constructible_from<
              BayesRecipe<mode_ad, PooledGaussianMethod>,
              SpikeSlabAD,
              VarianceBudget>);

// A configurable method accepts its own genetic spec or defaults that spec.
static_assert(std::constructible_from<
              BayesRecipe<mode_ad, UnpooledSpikeSlabMethod>,
              SpikeSlabAD,
              VarianceBudget>);
static_assert(!std::constructible_from<
              BayesRecipe<mode_ad, UnpooledSpikeSlabMethod>,
              ScaledMixtureAD,
              VarianceBudget>);
static_assert(std::constructible_from<
              BayesRecipe<mode_ad, UnpooledSpikeSlabMethod>,
              VarianceBudget>);
static_assert(std::constructible_from<
              BayesRecipe<mode_ad, DefaultScaledMixtureMethod>,
              ScaledMixtureAD,
              VarianceBudget>);
static_assert(std::constructible_from<
              BayesRecipe<mode_ad, DefaultJointSpikeSlabMethod>,
              JointSpikeSlabAD,
              VarianceBudget>);

// CD is a joint allocation across A and D, so it has no single-mode form.
static_assert(RecipeExists<mode_ad, DefaultJointSpikeSlabMethod>);
static_assert(!RecipeExists<mode_a, DefaultJointSpikeSlabMethod>);

static_assert(!std::same_as<
              BayesRecipe<mode_a, FixedPooledSpikeSlabMethod>,
              BayesRecipe<mode_a, PooledSpikeSlabMethod>>);
static_assert(
    std::same_as<
        std::remove_cvref_t<
            decltype(BayesRecipe<mode_a, FixedPooledSpikeSlabMethod>::defaults()
                         .genetic_spec())>,
        std::remove_cvref_t<
            decltype(BayesRecipe<mode_a, PooledSpikeSlabMethod>::defaults()
                         .genetic_spec())>>);

auto message_of(auto&& construct) -> std::string
{
    try
    {
        construct();
    }
    catch (const GelexException& error)
    {
        return error.what();
    }
    return {};
}

}  // namespace

TEST_CASE("BayesRecipe accepts a well-formed input", "[bayes][recipe]")
{
    const auto defaults = ScaledMixture{};
    const auto recipe = BayesRecipe<mode_ad, DefaultScaledMixtureMethod>{
        ScaledMixtureAD{ScaledMixture{}, ScaledMixture{}},
        VarianceBudget{{.additive = 0.4, .dominance = 0.05, .random = 0.05}},
    };

    REQUIRE(recipe.variance().share(GeneticMode::A) == 0.4);
    REQUIRE(recipe.variance().random() == 0.05);
    REQUIRE(
        recipe.genetic_spec().get<GeneticMode::D>().scales()
        == defaults.scales());
}

TEST_CASE(
    "BayesRecipe accepts default genetic specs with an explicit variance "
    "budget",
    "[bayes][recipe]")
{
    const auto recipe = BayesRecipe<mode_ad, DefaultScaledMixtureMethod>{
        VarianceBudget{{.additive = 0.4, .dominance = 0.05}}};
    const auto defaults = ScaledMixture{};

    REQUIRE(
        recipe.genetic_spec().get<GeneticMode::A>().probabilities()
        == defaults.probabilities());
    REQUIRE(
        recipe.genetic_spec().get<GeneticMode::D>().scales()
        == defaults.scales());
}

TEST_CASE(
    "BayesRecipe cross-checks its variance budget against its modes",
    "[bayes][recipe]")
{
    SECTION("a present mode needs a positive share")
    {
        const auto message = message_of(
            []
            {
                return BayesRecipe<mode_ad, PooledGaussianMethod>{
                    VarianceBudget{{.additive = 0.5}}};
            });

        REQUIRE_THAT(
            message,
            ContainsSubstring(
                "D variance share must be positive when the mode is present"));
    }

    SECTION("an absent mode needs a zero share")
    {
        const auto message = message_of(
            []
            {
                return BayesRecipe<mode_a, PooledGaussianMethod>{
                    VarianceBudget{{.additive = 0.5, .dominance = 0.2}}};
            });

        REQUIRE_THAT(
            message, ContainsSubstring("D variance share must be zero"));
    }
}

TEST_CASE("BayesRecipe::defaults fills in the mode defaults", "[bayes][recipe]")
{
    const auto additive_only
        = BayesRecipe<mode_a, PooledGaussianMethod>::defaults();

    REQUIRE(additive_only.variance().share(GeneticMode::A) == 0.5);
    REQUIRE(additive_only.variance().share(GeneticMode::D) == 0.0);
    REQUIRE(additive_only.variance().random() == 0.0);

    const auto both
        = BayesRecipe<mode_ad, DefaultScaledMixtureMethod>::defaults();
    const auto defaults = ScaledMixture{};

    REQUIRE(both.variance().share(GeneticMode::D) == 0.2);
    REQUIRE(
        both.genetic_spec().get<GeneticMode::A>().probabilities()
        == defaults.probabilities());
}

TEST_CASE("BayesRecipe::defaults covers every joint default", "[bayes][recipe]")
{
    const auto recipe
        = BayesRecipe<mode_ad, DefaultJointSpikeSlabMethod>::defaults();
    const auto defaults = JointSpikeSlab{};

    REQUIRE(
        recipe.genetic_spec().joint().probabilities()
        == defaults.probabilities());
    REQUIRE(
        recipe.genetic_spec()
            .mode_values()
            .get<GeneticMode::D>()
            .positive_probability()
        == 0.5);

    const auto magnitude_recipe
        = BayesRecipe<mode_ad, MagnitudeJointSpikeSlabMethod>::defaults();
    REQUIRE(
        magnitude_recipe.genetic_spec().joint().probabilities()
        == defaults.probabilities());
}
