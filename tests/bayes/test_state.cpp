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
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "gelex/bayes/detail/state_factory.h"
#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/joint_spike_slab.h"
#include "gelex/bayes/genetic/scaled_mixture.h"
#include "gelex/bayes/genetic/spike_slab.h"
#include "gelex/bayes/genetic_policy.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/variance/budget.h"
#include "gelex/data/fixed_design.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

#include "compact_genotype_fixture.h"
#include "random_design_fixture.h"

using Catch::Approx;
using gelex::BayesModel;
using gelex::BayesPrior;
using gelex::BayesRecipe;
using gelex::BayesState;
using gelex::GaussianPrior;
using gelex::GaussianSpec;
using gelex::GaussianState;
using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::HalfNormalPrior;
using gelex::HalfNormalSpec;
using gelex::HalfNormalState;
using gelex::JointModeValues;
using gelex::JointSpikeSlabPrior;
using gelex::JointSpikeSlabSpec;
using gelex::JointSpikeSlabState;
using gelex::MixtureWeightUpdate;
using gelex::ModeValues;
using gelex::ScaledMixturePrior;
using gelex::ScaledMixtureSpec;
using gelex::ScaledMixtureState;
using gelex::SpikeSlabPrior;
using gelex::SpikeSlabSpec;
using gelex::SpikeSlabState;
using gelex::VarianceBudget;
using gelex::VarianceLayout;

namespace
{

constexpr auto mode_a = GeneticModeSet{GeneticMode::A};
constexpr auto mode_ad = GeneticMode::A | GeneticMode::D;

using FixedUnpooledSpikeSlabSpecAD = gelex::HomogeneousModeValues<
    mode_ad,
    SpikeSlabSpec<VarianceLayout::Unpooled, MixtureWeightUpdate::Disabled>>;
using ScaledMixtureSpecA
    = gelex::HomogeneousModeValues<mode_a, ScaledMixtureSpec<>>;
using FixedJointSpikeSlabSpecAD = gelex::JointModeValues<
    gelex::ModeValues<mode_ad, gelex::GaussianSpec<>, gelex::HalfNormalSpec>,
    JointSpikeSlabSpec<MixtureWeightUpdate::Disabled>>;

using PooledGaussianPriorAD = ModeValues<
    mode_ad,
    GaussianPrior<VarianceLayout::Pooled>,
    GaussianPrior<VarianceLayout::Pooled>>;
using UnpooledGaussianPriorA
    = ModeValues<mode_a, GaussianPrior<VarianceLayout::Unpooled>>;
using PooledSpikeSlabPriorAD = ModeValues<
    mode_ad,
    SpikeSlabPrior<VarianceLayout::Pooled>,
    SpikeSlabPrior<VarianceLayout::Pooled>>;
using UnpooledSpikeSlabPriorAD = ModeValues<
    mode_ad,
    SpikeSlabPrior<VarianceLayout::Unpooled>,
    SpikeSlabPrior<VarianceLayout::Unpooled>>;
using FixedUnpooledSpikeSlabPriorAD = ModeValues<
    mode_ad,
    SpikeSlabPrior<VarianceLayout::Unpooled, MixtureWeightUpdate::Disabled>,
    SpikeSlabPrior<VarianceLayout::Unpooled, MixtureWeightUpdate::Disabled>>;
using HeterogeneousPriorAD = ModeValues<
    mode_ad,
    GaussianPrior<VarianceLayout::Pooled>,
    SpikeSlabPrior<VarianceLayout::Unpooled, MixtureWeightUpdate::Disabled>>;
using ScaledMixturePriorAD
    = ModeValues<mode_ad, ScaledMixturePrior<>, ScaledMixturePrior<>>;
using JointPrior = JointModeValues<
    ModeValues<mode_ad, GaussianPrior<VarianceLayout::Pooled>, HalfNormalPrior>,
    JointSpikeSlabPrior<>>;
using JointModeSpecs = ModeValues<mode_ad, GaussianSpec<>, HalfNormalSpec>;

static_assert(ScaledMixtureState::class_count == 5);
static_assert(JointSpikeSlabState::class_count == 4);
static_assert(ScaledMixtureState::component_count == 4);
static_assert(JointSpikeSlabState::component_count == 4);
static_assert(
    std::same_as<
        decltype(std::declval<const SpikeSlabState<VarianceLayout::Pooled>&>()
                     .assignments()),
        const Eigen::VectorX<std::uint8_t>&>);
static_assert(std::same_as<
              decltype(std::declval<const ScaledMixtureState&>().assignments()),
              const Eigen::VectorX<std::uint8_t>&>);
static_assert(
    std::same_as<
        decltype(std::declval<const JointSpikeSlabState&>().assignments()),
        const Eigen::VectorX<std::uint8_t>&>);

static_assert(std::same_as<
              gelex::detail::genetic_state_t<PooledGaussianPriorAD>,
              ModeValues<
                  mode_ad,
                  GaussianState<VarianceLayout::Pooled>,
                  GaussianState<VarianceLayout::Pooled>>>);
static_assert(std::same_as<
              gelex::detail::genetic_state_t<UnpooledGaussianPriorA>,
              ModeValues<mode_a, GaussianState<VarianceLayout::Unpooled>>>);
static_assert(std::same_as<
              gelex::detail::genetic_state_t<PooledSpikeSlabPriorAD>,
              ModeValues<
                  mode_ad,
                  SpikeSlabState<VarianceLayout::Pooled>,
                  SpikeSlabState<VarianceLayout::Pooled>>>);
static_assert(std::same_as<
              gelex::detail::genetic_state_t<UnpooledSpikeSlabPriorAD>,
              ModeValues<
                  mode_ad,
                  SpikeSlabState<VarianceLayout::Unpooled>,
                  SpikeSlabState<VarianceLayout::Unpooled>>>);
static_assert(std::same_as<
              gelex::detail::genetic_state_t<UnpooledSpikeSlabPriorAD>,
              gelex::detail::genetic_state_t<FixedUnpooledSpikeSlabPriorAD>>);
static_assert(std::same_as<
              gelex::detail::genetic_state_t<HeterogeneousPriorAD>,
              ModeValues<
                  mode_ad,
                  GaussianState<VarianceLayout::Pooled>,
                  SpikeSlabState<VarianceLayout::Unpooled>>>);
static_assert(std::same_as<
              gelex::detail::genetic_state_t<ScaledMixturePriorAD>,
              ModeValues<mode_ad, ScaledMixtureState, ScaledMixtureState>>);
static_assert(std::same_as<
              gelex::detail::genetic_state_t<JointPrior>,
              JointModeValues<
                  ModeValues<
                      mode_ad,
                      GaussianState<VarianceLayout::Pooled>,
                      HalfNormalState>,
                  JointSpikeSlabState>>);

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
    random.push_back(
        gelex::test::make_random_design(
            "batch",
            std::vector<std::string>{"batch_1", "batch_2"},
            Eigen::MatrixXd{{1.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}}));
    return BayesModel{
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        gelex::FixedDesign::make(3),
        std::move(random),
        std::move(genetic)};
}

}  // namespace

TEST_CASE(
    "genetic state preserves independent topology and initializes mode storage",
    "[bayes][state]")
{
    const auto model = make_model(mode_ad);
    const auto prior = gelex::make_prior(
        BayesRecipe<mode_ad, GaussianSpec<VarianceLayout::Pooled>>::defaults(),
        model);

    const auto state
        = gelex::detail::make_state(prior.genetic(), model.genetic());

    state.for_each(
        [&]<GeneticMode Mode>(const auto& mode_state)
        {
            STATIC_REQUIRE(mode_ad.contains(Mode));
            REQUIRE(mode_state.coefficients().size() == model.genetic().cols());
            REQUIRE(mode_state.coefficients().isZero());
            REQUIRE(
                mode_state.fitted_values().size() == model.genetic().rows());
            REQUIRE(mode_state.fitted_values().isZero());
        });
    REQUIRE(
        state.get<GeneticMode::A>().variance()
        == Approx(prior.genetic().get<GeneticMode::A>().variance.initial));
    REQUIRE(
        state.get<GeneticMode::D>().variance()
        == Approx(prior.genetic().get<GeneticMode::D>().variance.initial));
}

TEST_CASE(
    "unpooled Gaussian state expands the calibrated variance per marker",
    "[bayes][state]")
{
    const auto model = make_model(mode_a);
    const auto prior = gelex::make_prior(
        BayesRecipe<mode_a, GaussianSpec<VarianceLayout::Unpooled>>::defaults(),
        model);

    const auto state
        = gelex::detail::make_state(prior.genetic(), model.genetic());
    const auto& variance = state.get<GeneticMode::A>().variance();

    REQUIRE(variance.isApprox(
        Eigen::VectorXd::Constant(
            model.genetic().cols(),
            prior.genetic().get<GeneticMode::A>().variance.initial)));
}

TEST_CASE(
    "fixed spike-slab state initializes every mode probability",
    "[bayes][state]")
{
    const auto model = make_model(mode_ad);
    const auto recipe = BayesRecipe<mode_ad, FixedUnpooledSpikeSlabSpecAD>{
        FixedUnpooledSpikeSlabSpecAD{
            SpikeSlabSpec<
                VarianceLayout::Unpooled,
                MixtureWeightUpdate::Disabled>{0.05},
            SpikeSlabSpec<
                VarianceLayout::Unpooled,
                MixtureWeightUpdate::Disabled>{0.2},
        },
        VarianceBudget{{.additive = 0.4, .dominance = 0.1}}};
    const auto prior = gelex::make_prior(recipe, model);

    const auto state
        = gelex::detail::make_state(prior.genetic(), model.genetic());
    const auto& additive = state.get<GeneticMode::A>();
    const auto& dominance = state.get<GeneticMode::D>();

    REQUIRE(additive.probability() == 0.05);
    REQUIRE(dominance.probability() == 0.2);
    REQUIRE(additive.assignments().size() == model.genetic().cols());
    REQUIRE(additive.assignments().isZero());
    REQUIRE(dominance.assignments().isZero());
}

TEST_CASE(
    "scaled-mixture state owns assignment and component fitted caches",
    "[bayes][state]")
{
    const auto model = make_model(mode_a);
    const auto prior = gelex::make_prior(
        BayesRecipe<mode_a, ScaledMixtureSpecA>::defaults(), model);

    const auto state
        = gelex::detail::make_state(prior.genetic(), model.genetic());
    const auto& mode_state = state.get<GeneticMode::A>();
    const auto& family_state = mode_state;

    REQUIRE(family_state.assignments().size() == model.genetic().cols());
    REQUIRE(family_state.assignments().isZero());
    REQUIRE(
        family_state.probabilities()
        == prior.genetic().get<GeneticMode::A>().probabilities.initial);
    REQUIRE(family_state.fitted_values().rows() == model.genetic().rows());
    REQUIRE(
        family_state.fitted_values().cols()
        == static_cast<Eigen::Index>(ScaledMixtureState::component_count));
    REQUIRE(family_state.fitted_values().isZero());
}

TEST_CASE(
    "joint state keeps shared latent variables in the joint value",
    "[bayes][state]")
{
    const auto model = make_model(mode_ad);
    const auto recipe = BayesRecipe<mode_ad, FixedJointSpikeSlabSpecAD>{
        FixedJointSpikeSlabSpecAD{
            JointModeSpecs{GaussianSpec<>{}, HalfNormalSpec{}},
            JointSpikeSlabSpec<MixtureWeightUpdate::Disabled>{
                {0.8, 0.1, 0.05, 0.05}}},
        VarianceBudget{{.additive = 0.4, .dominance = 0.1}}};
    const auto prior = gelex::make_prior(recipe, model);

    const auto state
        = gelex::detail::make_state(prior.genetic(), model.genetic());
    const auto& joint = state.joint();

    REQUIRE(
        joint.probabilities() == recipe.genetic_spec().joint().probabilities());
    REQUIRE(joint.assignments().size() == model.genetic().cols());
    REQUIRE(joint.assignments().isZero());
    REQUIRE(joint.fitted_values().rows() == model.genetic().rows());
    REQUIRE(
        joint.fitted_values().cols()
        == static_cast<Eigen::Index>(JointSpikeSlabState::component_count));
    REQUIRE(joint.fitted_values().isZero());
    // state.mode_values().for_each(
    //     [&]<GeneticMode Mode>(const auto& mode_state)
    //     {
    //         STATIC_REQUIRE(mode_ad.contains(Mode));
    //         REQUIRE(
    //             mode_state.family_state.
    //             == model.genetic().rows());
    //         REQUIRE(mode_state.family_state.fitted_values().isZero());
    //     });
    const auto& dominance = state.mode_values().get<GeneticMode::D>();
    REQUIRE(dominance.probit_coefficients().isZero());
    REQUIRE(
        state.mode_values().get<GeneticMode::A>().variance()
        == Approx(prior.genetic()
                      .mode_values()
                      .get<GeneticMode::A>()
                      .variance.initial));
}

TEST_CASE(
    "genetic state rejects a design missing a required mode",
    "[bayes][state]")
{
    const auto full_model = make_model(mode_ad);
    const auto additive_model = make_model(mode_a);
    const auto prior = gelex::make_prior(
        BayesRecipe<mode_ad, GaussianSpec<VarianceLayout::Pooled>>::defaults(),
        full_model);

    REQUIRE_THROWS_AS(
        gelex::detail::make_state(prior.genetic(), additive_model.genetic()),
        gelex::GelexException);
}

TEST_CASE(
    "aggregate state initializes every mutable axis from model and prior",
    "[bayes][state]")
{
    const auto model = make_model_with_random();
    const auto recipe
        = BayesRecipe<mode_a, GaussianSpec<VarianceLayout::Pooled>>{
            VarianceBudget{{.additive = 0.4, .random = 0.1}}};
    const auto prior = gelex::make_prior(recipe, model);

    const auto state = gelex::make_state(prior, model);

    REQUIRE(state.fixed().coefficients.size() == model.fixed().X().cols());
    REQUIRE(state.fixed().coefficients.isZero());
    REQUIRE(state.random().size() == 1);
    REQUIRE(
        state.random().front().coefficients.size()
        == model.random().front().X().cols());
    REQUIRE(state.random().front().coefficients.isZero());
    REQUIRE(
        state.random().front().variance
        == Approx(prior.random().front().initial));
    REQUIRE(state.residual().adjusted_response.isApprox(model.phenotype()));
    REQUIRE(state.residual().variance == Approx(prior.residual().initial));

    const auto& genetic = state.genetic().get<GeneticMode::A>();
    REQUIRE(genetic.coefficients().isZero());
    REQUIRE(genetic.fitted_values().isZero());
}

TEST_CASE(
    "aggregate state rejects a prior made for another random topology",
    "[bayes][state]")
{
    const auto random_model = make_model_with_random();
    const auto prior = gelex::make_prior(
        BayesRecipe<mode_a, GaussianSpec<VarianceLayout::Pooled>>{
            VarianceBudget{{.additive = 0.4, .random = 0.1}}},
        random_model);
    const auto model_without_random = make_model(mode_a);

    REQUIRE_THROWS_AS(
        gelex::make_state(prior, model_without_random), gelex::GelexException);
}
