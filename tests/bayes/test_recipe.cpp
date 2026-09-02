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
#include "gelex/bayes/genetic/independent_topology.h"
#include "gelex/bayes/genetic/joint_topology.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/semantic_method.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/variance_budget.h"
#include "gelex/exception.h"
#include "gelex/types/genetic_mode.h"

using Catch::Matchers::ContainsSubstring;
using gelex::BayesMethod;
using gelex::BayesRecipe;
using gelex::BuiltinBayesRecipe;
using gelex::Gaussian;
using gelex::GaussianMethod;
using gelex::GelexException;
using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::IndependentTopology;
using gelex::JointSpikeSlab;
using gelex::JointSpikeSlabMethod;
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

using SpikeSlabA = IndependentTopology<mode_a, SpikeSlab>;
using SpikeSlabAD = IndependentTopology<mode_ad, SpikeSlab, SpikeSlab>;
using ScaledMixtureAD
    = IndependentTopology<mode_ad, ScaledMixture, ScaledMixture>;
using JointSpikeSlabAD = gelex::JointTopology<
    IndependentTopology<
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

template <typename Method, GeneticModeSet Modes>
concept RecipeExists = requires { typename BayesRecipe<Method, Modes>; };

template <typename Recipe>
concept HasMethodToken = requires { Recipe::method; };

static_assert(BayesRecipe<DefaultScaledMixtureMethod, mode_a>::modes == mode_a);
static_assert(!HasMethodToken<BayesRecipe<DefaultScaledMixtureMethod, mode_a>>);
static_assert(std::same_as<
              BayesRecipe<PooledGaussianMethod, mode_a>::method_type,
              PooledGaussianMethod>);
static_assert(std::same_as<
              BuiltinBayesRecipe<BayesMethod::RR, mode_a>,
              BayesRecipe<PooledGaussianMethod, mode_a>>);
static_assert(std::same_as<
              BuiltinBayesRecipe<BayesMethod::A, mode_a>,
              BayesRecipe<UnpooledGaussianMethod, mode_a>>);
static_assert(std::same_as<
              BuiltinBayesRecipe<BayesMethod::B, mode_a>,
              BayesRecipe<UnpooledSpikeSlabMethod, mode_a>>);
static_assert(std::same_as<
              BuiltinBayesRecipe<BayesMethod::C, mode_a>,
              BayesRecipe<PooledSpikeSlabMethod, mode_a>>);
static_assert(std::same_as<
              BuiltinBayesRecipe<BayesMethod::R, mode_a>,
              BayesRecipe<DefaultScaledMixtureMethod, mode_a>>);
static_assert(std::same_as<
              BuiltinBayesRecipe<BayesMethod::CD, mode_ad>,
              BayesRecipe<DefaultJointSpikeSlabMethod, mode_ad>>);

// B and C take identical parameters but stay distinct types: the shared and
// per-marker variance difference is intrinsic to the method.
static_assert(!std::same_as<
              BayesRecipe<UnpooledSpikeSlabMethod, mode_ad>,
              BayesRecipe<PooledSpikeSlabMethod, mode_ad>>);

// RR and A have no structural parameters: their prior shape is fixed by the
// data alone. Gaussian says so at the type level and costs no storage.
static_assert(
    std::same_as<
        std::remove_cvref_t<
            decltype(BayesRecipe<PooledGaussianMethod, mode_ad>::defaults()
                         .parameters())>,
        Gaussian>);
static_assert(
    std::same_as<
        std::remove_cvref_t<
            decltype(BayesRecipe<UnpooledGaussianMethod, mode_ad>::defaults()
                         .parameters())>,
        Gaussian>);
static_assert(
    sizeof(BayesRecipe<PooledGaussianMethod, mode_ad>)
    == sizeof(VarianceBudget));
static_assert(std::constructible_from<
              BayesRecipe<PooledGaussianMethod, mode_ad>,
              VarianceBudget>);
static_assert(!std::constructible_from<
              BayesRecipe<PooledGaussianMethod, mode_ad>,
              SpikeSlabAD,
              VarianceBudget>);

// A parameterized method takes its own shape and only its own shape.
static_assert(std::constructible_from<
              BayesRecipe<UnpooledSpikeSlabMethod, mode_ad>,
              SpikeSlabAD,
              VarianceBudget>);
static_assert(!std::constructible_from<
              BayesRecipe<UnpooledSpikeSlabMethod, mode_ad>,
              ScaledMixtureAD,
              VarianceBudget>);
static_assert(!std::constructible_from<
              BayesRecipe<UnpooledSpikeSlabMethod, mode_ad>,
              VarianceBudget>);
static_assert(std::constructible_from<
              BayesRecipe<DefaultScaledMixtureMethod, mode_ad>,
              ScaledMixtureAD,
              VarianceBudget>);
static_assert(std::constructible_from<
              BayesRecipe<DefaultJointSpikeSlabMethod, mode_ad>,
              JointSpikeSlabAD,
              VarianceBudget>);

// CD is a joint allocation across A and D, so it has no single-mode form.
static_assert(RecipeExists<DefaultJointSpikeSlabMethod, mode_ad>);
static_assert(!RecipeExists<DefaultJointSpikeSlabMethod, mode_a>);

static_assert(!std::same_as<
              BayesRecipe<FixedPooledSpikeSlabMethod, mode_a>,
              BayesRecipe<PooledSpikeSlabMethod, mode_a>>);
static_assert(
    std::same_as<
        std::remove_cvref_t<
            decltype(BayesRecipe<FixedPooledSpikeSlabMethod, mode_a>::defaults()
                         .parameters())>,
        std::remove_cvref_t<
            decltype(BayesRecipe<PooledSpikeSlabMethod, mode_a>::defaults()
                         .parameters())>>);

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
    const auto recipe = BayesRecipe<DefaultScaledMixtureMethod, mode_ad>{
        ScaledMixtureAD{ScaledMixture{}, ScaledMixture{}},
        VarianceBudget{{.additive = 0.4, .dominance = 0.05, .random = 0.05}},
    };

    REQUIRE(recipe.variance().share(GeneticMode::A) == 0.4);
    REQUIRE(recipe.variance().random() == 0.05);
    REQUIRE(
        recipe.parameters().get<GeneticMode::D>().scales()
        == defaults.scales());
}

TEST_CASE("BayesRecipe validates its variance budget", "[bayes][recipe]")
{
    const auto message = message_of(
        []
        {
            return BayesRecipe<PooledGaussianMethod, mode_a>{
                VarianceBudget{{.additive = 0.5, .dominance = 0.2}}};
        });

    REQUIRE_THAT(message, ContainsSubstring("D variance share must be zero"));
}

TEST_CASE("BayesRecipe::defaults fills in the mode defaults", "[bayes][recipe]")
{
    const auto additive_only
        = BayesRecipe<PooledGaussianMethod, mode_a>::defaults();

    REQUIRE(additive_only.variance().share(GeneticMode::A) == 0.5);
    REQUIRE(additive_only.variance().share(GeneticMode::D) == 0.0);
    REQUIRE(additive_only.variance().random() == 0.0);

    const auto both
        = BayesRecipe<DefaultScaledMixtureMethod, mode_ad>::defaults();
    const auto defaults = ScaledMixture{};

    REQUIRE(both.variance().share(GeneticMode::D) == 0.2);
    REQUIRE(
        both.parameters().get<GeneticMode::A>().probabilities()
        == defaults.probabilities());
}

TEST_CASE("BayesRecipe::defaults covers every joint default", "[bayes][recipe]")
{
    const auto recipe
        = BayesRecipe<DefaultJointSpikeSlabMethod, mode_ad>::defaults();
    const auto defaults = JointSpikeSlab{};

    REQUIRE(
        recipe.parameters().joint().probabilities()
        == defaults.probabilities());
    REQUIRE(
        recipe.parameters()
            .mode_values()
            .get<GeneticMode::D>()
            .positive_probability()
        == 0.5);

    const auto magnitude_recipe
        = BayesRecipe<MagnitudeJointSpikeSlabMethod, mode_ad>::defaults();
    REQUIRE(
        magnitude_recipe.parameters().joint().probabilities()
        == defaults.probabilities());
}
