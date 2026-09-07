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
#include "gelex/bayes/genetic/types.h"
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
using gelex::GaussianSpec;
using gelex::GelexException;
using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::JointSpikeSlabSpec;
using gelex::MixtureWeightUpdate;
using gelex::ModeValues;
using gelex::ScaledMixtureSpec;
using gelex::SpikeSlabSpec;
using gelex::VarianceBudget;
using gelex::VarianceLayout;

namespace
{

constexpr auto mode_a = GeneticModeSet{GeneticMode::A};
constexpr auto mode_ad = GeneticMode::A | GeneticMode::D;

using UnpooledSpikeSlabSpecAD = gelex::
    HomogeneousModeValues<mode_ad, SpikeSlabSpec<VarianceLayout::Unpooled>>;
using FixedPooledSpikeSlabSpecAD = gelex::HomogeneousModeValues<
    mode_ad,
    SpikeSlabSpec<VarianceLayout::Pooled, MixtureWeightUpdate::Disabled>>;
using ScaledMixtureSpecAD
    = gelex::HomogeneousModeValues<mode_ad, ScaledMixtureSpec<>>;
using JointSpikeSlabSpecAD = gelex::JointModeValues<
    gelex::ModeValues<mode_ad, gelex::GaussianSpec<>, gelex::HalfNormalSpec>,
    JointSpikeSlabSpec<>>;

using SpikeSlabAD = ModeValues<mode_ad, SpikeSlabSpec<>, SpikeSlabSpec<>>;
using ScaledMixtureAD
    = ModeValues<mode_ad, ScaledMixtureSpec<>, ScaledMixtureSpec<>>;
using JointSpikeSlabAD = gelex::JointModeValues<
    ModeValues<mode_ad, GaussianSpec<>, gelex::HalfNormalSpec>,
    JointSpikeSlabSpec<>>;

using PooledGaussian = GaussianSpec<VarianceLayout::Pooled>;
using UnpooledGaussian = GaussianSpec<VarianceLayout::Unpooled>;
using UnpooledSpikeSlabAD = UnpooledSpikeSlabSpecAD;
using FixedSpikeSlabAD = FixedPooledSpikeSlabSpecAD;

static_assert(std::same_as<
              BuiltinBayesRecipe<mode_a, BayesMethod::RR>,
              BayesRecipe<
                  mode_a,
                  gelex::HomogeneousModeValues<mode_a, PooledGaussian>>>);
static_assert(std::same_as<
              BuiltinBayesRecipe<mode_a, BayesMethod::A>,
              BayesRecipe<
                  mode_a,
                  gelex::HomogeneousModeValues<mode_a, UnpooledGaussian>>>);
static_assert(std::same_as<
              BuiltinBayesRecipe<mode_ad, BayesMethod::B>,
              BayesRecipe<mode_ad, UnpooledSpikeSlabAD>>);
static_assert(std::same_as<
              BuiltinBayesRecipe<mode_ad, BayesMethod::C>,
              BayesRecipe<mode_ad, SpikeSlabAD>>);
static_assert(std::same_as<
              BuiltinBayesRecipe<mode_ad, BayesMethod::R>,
              BayesRecipe<mode_ad, ScaledMixtureAD>>);
static_assert(std::same_as<
              BuiltinBayesRecipe<mode_ad, BayesMethod::CD>,
              BayesRecipe<mode_ad, JointSpikeSlabAD>>);
static_assert(BayesRecipe<mode_ad, ScaledMixtureAD>::modes == mode_ad);
static_assert(std::same_as<
              BayesRecipe<
                  mode_ad,
                  gelex::HomogeneousModeValues<mode_ad, UnpooledGaussian>>::
                  genetic_spec_type,
              gelex::HomogeneousModeValues<mode_ad, UnpooledGaussian>>);
static_assert(std::same_as<
              decltype(BayesRecipe<mode_ad, FixedSpikeSlabAD>::defaults()
                           .genetic_spec()),
              const FixedSpikeSlabAD&>);
static_assert(std::constructible_from<
              BayesRecipe<
                  mode_ad,
                  gelex::HomogeneousModeValues<mode_ad, PooledGaussian>>,
              VarianceBudget>);
static_assert(std::constructible_from<
              BayesRecipe<mode_ad, UnpooledSpikeSlabAD>,
              UnpooledSpikeSlabAD,
              VarianceBudget>);
static_assert(!std::constructible_from<
              BayesRecipe<mode_ad, UnpooledSpikeSlabAD>,
              SpikeSlabAD,
              VarianceBudget>);
static_assert(!std::constructible_from<
              BayesRecipe<mode_ad, FixedSpikeSlabAD>,
              SpikeSlabAD,
              VarianceBudget>);
static_assert(!std::constructible_from<
              BayesRecipe<
                  mode_ad,
                  gelex::HomogeneousModeValues<mode_ad, PooledGaussian>>,
              gelex::HomogeneousModeValues<mode_ad, UnpooledGaussian>,
              VarianceBudget>);

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
    const auto defaults = ScaledMixtureSpec<>{};
    const auto recipe = BayesRecipe<mode_ad, ScaledMixtureSpecAD>{
        ScaledMixtureSpecAD{ScaledMixtureSpec<>{}, ScaledMixtureSpec<>{}},
        VarianceBudget{{.additive = 0.4, .dominance = 0.05, .random = 0.05}},
    };

    REQUIRE(recipe.variance().genetic(GeneticMode::A) == 0.4);
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
    const auto recipe = BayesRecipe<mode_ad, ScaledMixtureSpecAD>{
        VarianceBudget{{.additive = 0.4, .dominance = 0.05}}};
    const auto defaults = ScaledMixtureSpec<>{};

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
    SECTION("a present mode needs a positive proportion")
    {
        const auto message = message_of(
            []
            {
                return BayesRecipe<
                    mode_ad,
                    gelex::HomogeneousModeValues<
                        mode_ad,
                        GaussianSpec<VarianceLayout::Pooled>>>{
                    VarianceBudget{{.additive = 0.5}}};
            });

        REQUIRE_THAT(
            message,
            ContainsSubstring(
                "D variance proportion must be positive when the mode is "
                "present"));
    }

    SECTION("an absent mode needs a zero proportion")
    {
        const auto message = message_of(
            []
            {
                return BayesRecipe<
                    mode_a,
                    gelex::HomogeneousModeValues<
                        mode_a,
                        GaussianSpec<VarianceLayout::Pooled>>>{
                    VarianceBudget{{.additive = 0.5, .dominance = 0.2}}};
            });

        REQUIRE_THAT(
            message, ContainsSubstring("D variance proportion must be zero"));
    }
}

TEST_CASE("BayesRecipe::defaults fills in the mode defaults", "[bayes][recipe]")
{
    const auto additive_only = BayesRecipe<
        mode_a,
        gelex::HomogeneousModeValues<
            mode_a,
            GaussianSpec<VarianceLayout::Pooled>>>::defaults();

    REQUIRE(additive_only.variance().genetic(GeneticMode::A) == 0.5);
    REQUIRE(additive_only.variance().genetic(GeneticMode::D) == 0.0);
    REQUIRE(additive_only.variance().random() == 0.0);

    const auto both = BayesRecipe<mode_ad, ScaledMixtureSpecAD>::defaults();
    const auto defaults = ScaledMixtureSpec<>{};

    REQUIRE(both.variance().genetic(GeneticMode::D) == 0.2);
    REQUIRE(
        both.genetic_spec().get<GeneticMode::A>().probabilities()
        == defaults.probabilities());
}

TEST_CASE("BayesRecipe::defaults covers every joint default", "[bayes][recipe]")
{
    const auto recipe = BayesRecipe<mode_ad, JointSpikeSlabSpecAD>::defaults();
    const auto defaults = JointSpikeSlabSpec<>{};

    REQUIRE(
        recipe.genetic_spec().joint().probabilities()
        == defaults.probabilities());
}
