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

#include "gelex/bayes/genetic/independent_topology.h"
#include "gelex/bayes/method.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/variance_budget.h"
#include "gelex/exception.h"
#include "gelex/types/genetic_mode.h"

using Catch::Matchers::ContainsSubstring;
using gelex::BayesMethod;
using gelex::BayesRecipe;
using gelex::GelexException;
using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::IndependentTopology;
using gelex::JointSpikeSlab;
using gelex::NoParameters;
using gelex::ScaledMixture;
using gelex::SpikeSlab;
using gelex::UpdatePolicy;
using gelex::VarianceBudget;

namespace
{

constexpr auto MODE_A = GeneticModeSet{GeneticMode::A};
constexpr auto MODE_AD = GeneticMode::A | GeneticMode::D;

using SpikeSlabA = IndependentTopology<MODE_A, SpikeSlab>;
using SpikeSlabAD = IndependentTopology<MODE_AD, SpikeSlab>;
using ScaledMixtureAD = IndependentTopology<MODE_AD, ScaledMixture>;

const auto BUDGET_AD = VarianceBudget{{.additive = 0.5, .dominance = 0.2}};

template <BayesMethod Method, GeneticModeSet Modes>
concept RecipeExists = requires { typename BayesRecipe<Method, Modes>; };

// A recipe carries the method and mode set it was chosen with, so downstream
// compilation reads them off the type instead of a runtime tag.
static_assert(BayesRecipe<BayesMethod::R, MODE_A>::method == BayesMethod::R);
static_assert(BayesRecipe<BayesMethod::R, MODE_A>::modes == MODE_A);

// B and C take identical parameters but stay distinct types: the shared and
// per-marker variance difference is intrinsic to the method.
static_assert(!std::same_as<
              BayesRecipe<BayesMethod::B, MODE_AD>,
              BayesRecipe<BayesMethod::C, MODE_AD>>);

// RR and A have no structural parameters: their prior shape is fixed by the
// data alone. NoParameters says so at the type level and costs no storage.
static_assert(std::same_as<
              gelex::method_parameters_t<BayesMethod::RR, MODE_AD>,
              NoParameters>);
static_assert(std::same_as<
              gelex::method_parameters_t<BayesMethod::A, MODE_AD>,
              NoParameters>);
static_assert(
    sizeof(BayesRecipe<BayesMethod::RR, MODE_AD>) == sizeof(VarianceBudget));
static_assert(std::constructible_from<
              BayesRecipe<BayesMethod::RR, MODE_AD>,
              VarianceBudget>);
static_assert(!std::constructible_from<
              BayesRecipe<BayesMethod::RR, MODE_AD>,
              SpikeSlabAD,
              VarianceBudget>);

// A parameterized method takes its own shape and only its own shape.
static_assert(std::constructible_from<
              BayesRecipe<BayesMethod::B, MODE_AD>,
              SpikeSlabAD,
              VarianceBudget>);
static_assert(!std::constructible_from<
              BayesRecipe<BayesMethod::B, MODE_AD>,
              ScaledMixtureAD,
              VarianceBudget>);
static_assert(!std::constructible_from<
              BayesRecipe<BayesMethod::B, MODE_AD>,
              VarianceBudget>);
static_assert(std::constructible_from<
              BayesRecipe<BayesMethod::R, MODE_AD>,
              ScaledMixtureAD,
              VarianceBudget>);
static_assert(std::constructible_from<
              BayesRecipe<BayesMethod::CD, MODE_AD>,
              JointSpikeSlab,
              VarianceBudget>);

// CD is a joint allocation across A and D, so it has no single-mode form.
static_assert(RecipeExists<BayesMethod::CD, MODE_AD>);
static_assert(!RecipeExists<BayesMethod::CD, MODE_A>);
static_assert(!RecipeExists<BayesMethod::RR, GeneticModeSet{}>);

// Whether a parameter is sampled is a value, not part of the type.
static_assert(std::same_as<
              decltype(BayesRecipe<BayesMethod::C, MODE_A>{
                  SpikeSlabA{SpikeSlab{
                      .probability = {.update = UpdatePolicy::Fixed}}},
                  VarianceBudget{{.additive = 0.5}}}),
              decltype(BayesRecipe<BayesMethod::C, MODE_A>{
                  SpikeSlabA{SpikeSlab{}},
                  VarianceBudget{{.additive = 0.5}}})>);

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
    const auto recipe = BayesRecipe<BayesMethod::R, MODE_AD>{
        ScaledMixtureAD{ScaledMixture{}, ScaledMixture{}},
        VarianceBudget{{.additive = 0.4, .dominance = 0.05, .random = 0.05}},
    };

    REQUIRE(recipe.variance().share(GeneticMode::A) == 0.4);
    REQUIRE(recipe.variance().random() == 0.05);
    REQUIRE(
        recipe.parameters().get<GeneticMode::D>().scales
        == ScaledMixture{}.scales);
}

TEST_CASE("BayesRecipe validates its inputs on construction", "[bayes][recipe]")
{
    SECTION("the method parameters are checked")
    {
        const auto message = message_of(
            []
            {
                return BayesRecipe<BayesMethod::B, MODE_AD>{
                    SpikeSlabAD{
                        SpikeSlab{.probability = {.initial = 1.5}},
                        SpikeSlab{},
                    },
                    BUDGET_AD,
                };
            });

        REQUIRE_THAT(
            message, ContainsSubstring("A inclusion probability must lie"));
    }

    SECTION("the variance budget is checked against the recipe's own modes")
    {
        const auto message = message_of(
            []
            {
                return BayesRecipe<BayesMethod::RR, MODE_A>{
                    VarianceBudget{{.additive = 0.5, .dominance = 0.2}}};
            });

        REQUIRE_THAT(
            message, ContainsSubstring("D variance share must be zero"));
    }
}

TEST_CASE("BayesRecipe::defaults fills in the mode defaults", "[bayes][recipe]")
{
    const auto additive_only = BayesRecipe<BayesMethod::RR, MODE_A>::defaults();

    REQUIRE(additive_only.variance().share(GeneticMode::A) == 0.5);
    REQUIRE(additive_only.variance().share(GeneticMode::D) == 0.0);
    REQUIRE(additive_only.variance().random() == 0.0);

    const auto both = BayesRecipe<BayesMethod::R, MODE_AD>::defaults();

    REQUIRE(both.variance().share(GeneticMode::D) == 0.2);
    REQUIRE(
        both.parameters().get<GeneticMode::A>().probabilities.initial
        == ScaledMixture{}.probabilities.initial);
    REQUIRE(
        both.parameters().get<GeneticMode::D>().probabilities.update
        == UpdatePolicy::Sampled);
}

TEST_CASE("BayesRecipe::defaults covers every joint default", "[bayes][recipe]")
{
    const auto recipe = BayesRecipe<BayesMethod::CD, MODE_AD>::defaults();

    REQUIRE(
        recipe.parameters().positive_probability.initial
        == JointSpikeSlab{}.positive_probability.initial);
    REQUIRE(
        recipe.parameters().positive_probability.update
        == UpdatePolicy::Sampled);
}
