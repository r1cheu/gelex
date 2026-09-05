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

#include "gelex/bayes/builtin_method.h"
#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/joint_spike_slab.h"
#include "gelex/bayes/genetic/scaled_mixture.h"
#include "gelex/bayes/genetic/spike_slab.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/variance/budget.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

using Catch::Matchers::ContainsSubstring;
using gelex::BayesMethod;
using gelex::BayesRecipe;
using gelex::BuiltinBayesRecipe;
using gelex::Gaussian;
using gelex::GaussianFamily;
using gelex::GelexException;
using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::JointSpikeSlab;
using gelex::JointSpikeSlabFamily;
using gelex::MixtureWeightUpdate;
using gelex::ModeValues;
using gelex::ScaledMixture;
using gelex::ScaledMixtureFamily;
using gelex::SpikeSlab;
using gelex::SpikeSlabFamily;
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
    ModeValues<mode_ad, Gaussian, gelex::HalfNormal>,
    JointSpikeSlab>;
using PooledGaussianFamily = GaussianFamily<VarianceLayout::Pooled>;
using UnpooledGaussianFamily = GaussianFamily<VarianceLayout::Unpooled>;
using PooledSpikeSlabFamily = SpikeSlabFamily<VarianceLayout::Pooled>;
using UnpooledSpikeSlabFamily = SpikeSlabFamily<VarianceLayout::Unpooled>;
using FixedPooledSpikeSlabFamily
    = SpikeSlabFamily<VarianceLayout::Pooled, MixtureWeightUpdate::Disabled>;
using DefaultScaledMixtureFamily = ScaledMixtureFamily<>;
using DefaultJointSpikeSlabFamily = JointSpikeSlabFamily<>;

template <GeneticModeSet Modes, typename Family>
concept RecipeExists = requires { typename BayesRecipe<Modes, Family>; };

template <typename Recipe>
concept HasFamilyToken = requires { Recipe::family; };

static_assert(BayesRecipe<mode_a, DefaultScaledMixtureFamily>::modes == mode_a);
static_assert(!HasFamilyToken<BayesRecipe<mode_a, DefaultScaledMixtureFamily>>);
static_assert(std::same_as<
              BayesRecipe<mode_a, PooledGaussianFamily>::family_type,
              PooledGaussianFamily>);
static_assert(std::same_as<
              BuiltinBayesRecipe<mode_a, BayesMethod::RR>,
              BayesRecipe<mode_a, PooledGaussianFamily>>);
static_assert(std::same_as<
              BuiltinBayesRecipe<mode_a, BayesMethod::A>,
              BayesRecipe<mode_a, UnpooledGaussianFamily>>);
static_assert(std::same_as<
              BuiltinBayesRecipe<mode_a, BayesMethod::B>,
              BayesRecipe<mode_a, UnpooledSpikeSlabFamily>>);
static_assert(std::same_as<
              BuiltinBayesRecipe<mode_a, BayesMethod::C>,
              BayesRecipe<mode_a, PooledSpikeSlabFamily>>);
static_assert(std::same_as<
              BuiltinBayesRecipe<mode_a, BayesMethod::R>,
              BayesRecipe<mode_a, DefaultScaledMixtureFamily>>);
static_assert(std::same_as<
              BuiltinBayesRecipe<mode_ad, BayesMethod::CD>,
              BayesRecipe<mode_ad, DefaultJointSpikeSlabFamily>>);

// B and C take identical genetic specs but stay distinct types: the shared and
// per-marker variance difference is intrinsic to the family.
static_assert(!std::same_as<
              BayesRecipe<mode_ad, UnpooledSpikeSlabFamily>,
              BayesRecipe<mode_ad, PooledSpikeSlabFamily>>);

// RR and A have empty genetic specs: their prior shape is fixed by the data
// alone. Gaussian says so at the type level and costs no storage.
static_assert(
    std::same_as<
        typename BayesRecipe<mode_ad, PooledGaussianFamily>::genetic_spec_type,
        Gaussian>);
static_assert(
    std::same_as<
        std::remove_cvref_t<
            decltype(BayesRecipe<mode_ad, UnpooledGaussianFamily>::defaults()
                         .genetic_spec())>,
        Gaussian>);
static_assert(
    sizeof(BayesRecipe<mode_ad, PooledGaussianFamily>)
    == sizeof(VarianceBudget));
static_assert(std::constructible_from<
              BayesRecipe<mode_ad, PooledGaussianFamily>,
              VarianceBudget>);
static_assert(!std::constructible_from<
              BayesRecipe<mode_ad, PooledGaussianFamily>,
              SpikeSlabAD,
              VarianceBudget>);

// A configurable family accepts its own genetic spec or defaults that spec.
static_assert(std::constructible_from<
              BayesRecipe<mode_ad, UnpooledSpikeSlabFamily>,
              SpikeSlabAD,
              VarianceBudget>);
static_assert(!std::constructible_from<
              BayesRecipe<mode_ad, UnpooledSpikeSlabFamily>,
              ScaledMixtureAD,
              VarianceBudget>);
static_assert(std::constructible_from<
              BayesRecipe<mode_ad, UnpooledSpikeSlabFamily>,
              VarianceBudget>);
static_assert(std::constructible_from<
              BayesRecipe<mode_ad, DefaultScaledMixtureFamily>,
              ScaledMixtureAD,
              VarianceBudget>);
static_assert(std::constructible_from<
              BayesRecipe<mode_ad, DefaultJointSpikeSlabFamily>,
              JointSpikeSlabAD,
              VarianceBudget>);

// CD is a joint allocation across A and D, so it has no single-mode form.
static_assert(RecipeExists<mode_ad, DefaultJointSpikeSlabFamily>);
static_assert(!RecipeExists<mode_a, DefaultJointSpikeSlabFamily>);

static_assert(!std::same_as<
              BayesRecipe<mode_a, FixedPooledSpikeSlabFamily>,
              BayesRecipe<mode_a, PooledSpikeSlabFamily>>);
static_assert(
    std::same_as<
        std::remove_cvref_t<
            decltype(BayesRecipe<mode_a, FixedPooledSpikeSlabFamily>::defaults()
                         .genetic_spec())>,
        std::remove_cvref_t<
            decltype(BayesRecipe<mode_a, PooledSpikeSlabFamily>::defaults()
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
    const auto recipe = BayesRecipe<mode_ad, DefaultScaledMixtureFamily>{
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
    const auto recipe = BayesRecipe<mode_ad, DefaultScaledMixtureFamily>{
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
                return BayesRecipe<mode_ad, PooledGaussianFamily>{
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
                return BayesRecipe<mode_a, PooledGaussianFamily>{
                    VarianceBudget{{.additive = 0.5, .dominance = 0.2}}};
            });

        REQUIRE_THAT(
            message, ContainsSubstring("D variance share must be zero"));
    }
}

TEST_CASE("BayesRecipe::defaults fills in the mode defaults", "[bayes][recipe]")
{
    const auto additive_only
        = BayesRecipe<mode_a, PooledGaussianFamily>::defaults();

    REQUIRE(additive_only.variance().share(GeneticMode::A) == 0.5);
    REQUIRE(additive_only.variance().share(GeneticMode::D) == 0.0);
    REQUIRE(additive_only.variance().random() == 0.0);

    const auto both
        = BayesRecipe<mode_ad, DefaultScaledMixtureFamily>::defaults();
    const auto defaults = ScaledMixture{};

    REQUIRE(both.variance().share(GeneticMode::D) == 0.2);
    REQUIRE(
        both.genetic_spec().get<GeneticMode::A>().probabilities()
        == defaults.probabilities());
}

TEST_CASE("BayesRecipe::defaults covers every joint default", "[bayes][recipe]")
{
    const auto recipe
        = BayesRecipe<mode_ad, DefaultJointSpikeSlabFamily>::defaults();
    const auto defaults = JointSpikeSlab{};

    REQUIRE(
        recipe.genetic_spec().joint().probabilities()
        == defaults.probabilities());
}
