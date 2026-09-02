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
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <concepts>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "gelex/bayes/detail/genetic_state_compilation.h"
#include "gelex/bayes/genetic/component_layout.h"
#include "gelex/bayes/genetic/independent_topology.h"
#include "gelex/bayes/genetic/joint_topology.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior_compilation.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/semantic_method.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/variance_budget.h"
#include "gelex/exception.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_mode.h"

#include "compact_genotype_fixture.h"

using Catch::Approx;
using gelex::BayesModel;
using gelex::BayesPrior;
using gelex::BayesRecipe;
using gelex::BayesState;
using gelex::GaussianMethod;
using gelex::GaussianPrior;
using gelex::GaussianState;
using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::GeneticModeState;
using gelex::IndependentTopology;
using gelex::JointSpikeSlabMethod;
using gelex::JointSpikeSlabPrior;
using gelex::JointSpikeSlabState;
using gelex::JointTopology;
using gelex::JointZeroInflatedComponentLayout;
using gelex::ScaledMixture;
using gelex::ScaledMixtureMethod;
using gelex::ScaledMixturePrior;
using gelex::ScaledMixtureState;
using gelex::SingleComponentLayout;
using gelex::SpikeSlab;
using gelex::SpikeSlabMethod;
using gelex::SpikeSlabPrior;
using gelex::SpikeSlabState;
using gelex::UpdatePolicy;
using gelex::VarianceBudget;
using gelex::VarianceLayout;
using gelex::ZeroInflatedComponentLayout;

namespace
{

constexpr auto mode_a = GeneticModeSet{GeneticMode::A};
constexpr auto mode_ad = GeneticMode::A | GeneticMode::D;

using PooledGaussianPriorAD
    = IndependentTopology<mode_ad, GaussianPrior<VarianceLayout::Pooled>>;
using PooledGaussianPriorA
    = IndependentTopology<mode_a, GaussianPrior<VarianceLayout::Pooled>>;
using UnpooledGaussianPriorA
    = IndependentTopology<mode_a, GaussianPrior<VarianceLayout::Unpooled>>;
using PooledSpikeSlabPriorAD
    = IndependentTopology<mode_ad, SpikeSlabPrior<VarianceLayout::Pooled>>;
using UnpooledSpikeSlabPriorAD
    = IndependentTopology<mode_ad, SpikeSlabPrior<VarianceLayout::Unpooled>>;
using FixedUnpooledSpikeSlabPriorAD = IndependentTopology<
    mode_ad,
    SpikeSlabPrior<VarianceLayout::Unpooled, UpdatePolicy::Fixed>>;
using ScaledMixturePriorAD = IndependentTopology<mode_ad, ScaledMixturePrior<>>;
using JointPrior = JointTopology<
    GaussianPrior<VarianceLayout::Pooled>,
    JointSpikeSlabPrior<>>;

static_assert(std::same_as<
              decltype(SpikeSlabState<VarianceLayout::Pooled>::assignment),
              Eigen::VectorX<std::uint8_t>>);
static_assert(std::same_as<
              decltype(ScaledMixtureState::assignment),
              Eigen::VectorX<std::uint8_t>>);
static_assert(std::same_as<
              decltype(JointSpikeSlabState::assignment),
              Eigen::VectorX<std::uint8_t>>);
static_assert(std::same_as<
              decltype(JointSpikeSlabState::dominance_sign),
              Eigen::VectorX<std::uint8_t>>);

static_assert(std::same_as<
              gelex::detail::genetic_state_t<PooledGaussianPriorAD>,
              IndependentTopology<
                  mode_ad,
                  GeneticModeState<
                      GaussianState<VarianceLayout::Pooled>,
                      SingleComponentLayout>>>);
static_assert(std::same_as<
              gelex::detail::genetic_state_t<UnpooledGaussianPriorA>,
              IndependentTopology<
                  mode_a,
                  GeneticModeState<
                      GaussianState<VarianceLayout::Unpooled>,
                      SingleComponentLayout>>>);
static_assert(std::same_as<
              gelex::detail::genetic_state_t<PooledSpikeSlabPriorAD>,
              IndependentTopology<
                  mode_ad,
                  GeneticModeState<
                      SpikeSlabState<VarianceLayout::Pooled>,
                      ZeroInflatedComponentLayout<2>>>>);
static_assert(std::same_as<
              gelex::detail::genetic_state_t<UnpooledSpikeSlabPriorAD>,
              IndependentTopology<
                  mode_ad,
                  GeneticModeState<
                      SpikeSlabState<VarianceLayout::Unpooled>,
                      ZeroInflatedComponentLayout<2>>>>);
static_assert(std::same_as<
              gelex::detail::genetic_state_t<UnpooledSpikeSlabPriorAD>,
              gelex::detail::genetic_state_t<FixedUnpooledSpikeSlabPriorAD>>);
static_assert(std::same_as<
              gelex::detail::genetic_state_t<ScaledMixturePriorAD>,
              IndependentTopology<
                  mode_ad,
                  GeneticModeState<
                      ScaledMixtureState,
                      ZeroInflatedComponentLayout<
                          ScaledMixturePrior<>::class_count>>>>);
static_assert(std::same_as<
              gelex::detail::genetic_state_t<JointPrior>,
              JointTopology<
                  GeneticModeState<
                      GaussianState<VarianceLayout::Pooled>,
                      JointZeroInflatedComponentLayout>,
                  JointSpikeSlabState>>);
static_assert(std::same_as<
              gelex::bayes_state_t<BayesPrior<PooledGaussianPriorA>>,
              BayesState<PooledGaussianPriorA>>);

using PooledGaussianMethod = GaussianMethod<VarianceLayout::Pooled>;
using UnpooledGaussianMethod = GaussianMethod<VarianceLayout::Unpooled>;
using FixedUnpooledSpikeSlabMethod
    = SpikeSlabMethod<VarianceLayout::Unpooled, UpdatePolicy::Fixed>;
using DefaultScaledMixtureMethod = ScaledMixtureMethod<>;
using MixedJointSpikeSlabMethod
    = JointSpikeSlabMethod<UpdatePolicy::Fixed, UpdatePolicy::Sampled>;
using SpikeSlabAD = IndependentTopology<mode_ad, SpikeSlab>;

auto make_model(GeneticModeSet modes) -> BayesModel
{
    return gelex::test::make_compact_model(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}},
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        modes);
}

auto make_model_with_random() -> BayesModel
{
    auto genetic = gelex::test::make_genetic_design(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}}, mode_a);
    std::vector<gelex::bayes::RandomDesign> random;
    random.emplace_back(
        "batch",
        std::vector<std::string>{"batch_1", "batch_2"},
        Eigen::MatrixXd{{1.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}});
    return BayesModel{
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        gelex::FixedDesign::make(3),
        std::move(random),
        std::move(genetic)};
}

}  // namespace

TEST_CASE(
    "genetic state preserves independent topology and initializes mode storage",
    "[bayes][state_compilation]")
{
    const auto model = make_model(mode_ad);
    const auto prior = gelex::compile(
        BayesRecipe<PooledGaussianMethod, mode_ad>::defaults(), model);

    const auto state
        = gelex::detail::make_genetic_state(prior.genetic(), model.genetic());

    for (const auto& [mode, mode_state] : state.each())
    {
        REQUIRE(mode_ad.contains(mode));
        REQUIRE(mode_state.coefficients.size() == model.genetic().cols());
        REQUIRE(mode_state.coefficients.isZero());
        REQUIRE(
            mode_state.component_fitted_values.rows()
            == model.genetic().rows());
        REQUIRE(mode_state.component_fitted_values.cols() == 1);
        REQUIRE(mode_state.component_fitted_values.isZero());
    }
    REQUIRE(
        state.get<GeneticMode::A>().family_state.variance
        == Approx(
            prior.genetic().get<GeneticMode::A>().variance.initial_value()));
    REQUIRE(
        state.get<GeneticMode::D>().family_state.variance
        == Approx(
            prior.genetic().get<GeneticMode::D>().variance.initial_value()));
}

TEST_CASE(
    "unpooled Gaussian state expands the calibrated variance per marker",
    "[bayes][state_compilation]")
{
    const auto model = make_model(mode_a);
    const auto prior = gelex::compile(
        BayesRecipe<UnpooledGaussianMethod, mode_a>::defaults(), model);

    const auto state
        = gelex::detail::make_genetic_state(prior.genetic(), model.genetic());
    const auto& variance = state.get<GeneticMode::A>().family_state.variance;

    REQUIRE(variance.isApprox(
        Eigen::VectorXd::Constant(
            model.genetic().cols(),
            prior.genetic().get<GeneticMode::A>().variance.initial_value())));
}

TEST_CASE(
    "fixed spike-slab state initializes every mode probability",
    "[bayes][state_compilation]")
{
    const auto model = make_model(mode_ad);
    const auto recipe = BayesRecipe<FixedUnpooledSpikeSlabMethod, mode_ad>{
        SpikeSlabAD{
            SpikeSlab{.probability = 0.05},
            SpikeSlab{.probability = 0.2},
        },
        VarianceBudget{{.additive = 0.4, .dominance = 0.1}}};
    const auto prior = gelex::compile(recipe, model);

    const auto state
        = gelex::detail::make_genetic_state(prior.genetic(), model.genetic());
    const auto& additive = state.get<GeneticMode::A>().family_state;
    const auto& dominance = state.get<GeneticMode::D>().family_state;

    REQUIRE(additive.probability == 0.05);
    REQUIRE(dominance.probability == 0.2);
    REQUIRE(additive.assignment.size() == model.genetic().cols());
    REQUIRE(additive.assignment.isZero());
    REQUIRE(dominance.assignment.isZero());
}

TEST_CASE(
    "scaled-mixture state owns assignment and component fitted caches",
    "[bayes][state_compilation]")
{
    const auto model = make_model(mode_a);
    const auto prior = gelex::compile(
        BayesRecipe<DefaultScaledMixtureMethod, mode_a>::defaults(), model);

    const auto state
        = gelex::detail::make_genetic_state(prior.genetic(), model.genetic());
    const auto& mode_state = state.get<GeneticMode::A>();
    const auto& family_state = mode_state.family_state;

    REQUIRE(family_state.assignment.size() == model.genetic().cols());
    REQUIRE(family_state.assignment.isZero());
    REQUIRE(
        family_state.probabilities
        == prior.genetic().get<GeneticMode::A>().probabilities.initial);
    REQUIRE(
        mode_state.component_fitted_values.rows() == model.genetic().rows());
    REQUIRE(
        mode_state.component_fitted_values.cols()
        == static_cast<Eigen::Index>(
            ScaledMixturePrior<>::component_layout::component_count));
    REQUIRE(mode_state.component_fitted_values.isZero());
}

TEST_CASE(
    "joint state keeps shared latent variables in the joint value",
    "[bayes][state_compilation]")
{
    const auto model = make_model(mode_ad);
    const auto recipe = BayesRecipe<MixedJointSpikeSlabMethod, mode_ad>{
        gelex::JointSpikeSlab{
            .probabilities = {0.8, 0.1, 0.05, 0.05},
            .positive_probability = 0.6},
        VarianceBudget{{.additive = 0.4, .dominance = 0.1}}};
    const auto prior = gelex::compile(recipe, model);

    const auto state
        = gelex::detail::make_genetic_state(prior.genetic(), model.genetic());
    const auto& joint = state.joint();

    REQUIRE(joint.probabilities == recipe.parameters().probabilities);
    REQUIRE(joint.positive_probability == 0.6);
    REQUIRE(joint.assignment.size() == model.genetic().cols());
    REQUIRE(joint.assignment.isZero());
    REQUIRE(joint.dominance_sign.size() == model.genetic().cols());
    REQUIRE(joint.dominance_sign.isZero());
    for (const auto& [mode, mode_state] : state.mode_values().each())
    {
        REQUIRE(mode_ad.contains(mode));
        REQUIRE(
            mode_state.component_fitted_values.rows()
            == model.genetic().rows());
        REQUIRE(
            mode_state.component_fitted_values.cols()
            == static_cast<Eigen::Index>(
                JointZeroInflatedComponentLayout::component_count));
        REQUIRE(mode_state.component_fitted_values.isZero());
    }
    REQUIRE(
        state.mode_values().get<GeneticMode::A>().family_state.variance
        == Approx(prior.genetic()
                      .mode_values()
                      .get<GeneticMode::A>()
                      .variance.initial_value()));
}

TEST_CASE(
    "genetic state rejects a design missing a required mode",
    "[bayes][state_compilation]")
{
    const auto full_model = make_model(mode_ad);
    const auto additive_model = make_model(mode_a);
    const auto prior = gelex::compile(
        BayesRecipe<PooledGaussianMethod, mode_ad>::defaults(), full_model);

    REQUIRE_THROWS_AS(
        gelex::detail::make_genetic_state(
            prior.genetic(), additive_model.genetic()),
        gelex::GelexException);
}

TEST_CASE(
    "aggregate state initializes every mutable axis from model and prior",
    "[bayes][state_compilation]")
{
    const auto model = make_model_with_random();
    const auto recipe = BayesRecipe<PooledGaussianMethod, mode_a>{
        VarianceBudget{{.additive = 0.4, .random = 0.1}}};
    const auto prior = gelex::compile(recipe, model);

    const auto state = gelex::make_state(prior, model);

    REQUIRE(state.fixed().coefficients.size() == model.fixed().X.cols());
    REQUIRE(state.fixed().coefficients.isZero());
    REQUIRE(state.random().size() == 1);
    REQUIRE(
        state.random().front().coefficients.size()
        == model.random().front().X.cols());
    REQUIRE(state.random().front().coefficients.isZero());
    REQUIRE(
        state.random().front().variance
        == Approx(prior.random().front().initial_value()));
    REQUIRE(state.residual().adjusted_response.isApprox(model.phenotype()));
    REQUIRE(
        state.residual().variance == Approx(prior.residual().initial_value()));

    const auto& genetic = state.genetic().get<GeneticMode::A>();
    REQUIRE(genetic.coefficients.isZero());
    REQUIRE(genetic.component_fitted_values.isZero());
}

TEST_CASE(
    "aggregate state rejects a prior compiled for another random topology",
    "[bayes][state_compilation]")
{
    const auto random_model = make_model_with_random();
    const auto prior = gelex::compile(
        BayesRecipe<PooledGaussianMethod, mode_a>{
            VarianceBudget{{.additive = 0.4, .random = 0.1}}},
        random_model);
    const auto model_without_random = make_model(mode_a);

    REQUIRE_THROWS_AS(
        gelex::make_state(prior, model_without_random), gelex::GelexException);
}
