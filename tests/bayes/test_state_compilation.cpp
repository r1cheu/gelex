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

#include "gelex/bayes/detail/genetic_state_compilation.h"
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
using gelex::Gaussian;
using gelex::GaussianMethod;
using gelex::GaussianPrior;
using gelex::GaussianState;
using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::GeneticModeState;
using gelex::HalfNormal;
using gelex::HalfNormalAsymmetry;
using gelex::HalfNormalPrior;
using gelex::HalfNormalState;
using gelex::IndependentTopology;
using gelex::JointSpikeSlab;
using gelex::JointSpikeSlabMethod;
using gelex::JointSpikeSlabPrior;
using gelex::JointSpikeSlabState;
using gelex::JointTopology;
using gelex::ScaledMixture;
using gelex::ScaledMixtureMethod;
using gelex::ScaledMixturePrior;
using gelex::ScaledMixtureState;
using gelex::SpikeSlab;
using gelex::SpikeSlabMethod;
using gelex::SpikeSlabPrior;
using gelex::SpikeSlabState;
using gelex::UpdatePolicy;
using gelex::VarianceBudget;
using gelex::VarianceLayout;

namespace
{

constexpr auto mode_a = GeneticModeSet{GeneticMode::A};
constexpr auto mode_ad = GeneticMode::A | GeneticMode::D;

using PooledGaussianPriorAD = IndependentTopology<
    mode_ad,
    GaussianPrior<VarianceLayout::Pooled>,
    GaussianPrior<VarianceLayout::Pooled>>;
using PooledGaussianPriorA
    = IndependentTopology<mode_a, GaussianPrior<VarianceLayout::Pooled>>;
using UnpooledGaussianPriorA
    = IndependentTopology<mode_a, GaussianPrior<VarianceLayout::Unpooled>>;
using PooledSpikeSlabPriorAD = IndependentTopology<
    mode_ad,
    SpikeSlabPrior<VarianceLayout::Pooled>,
    SpikeSlabPrior<VarianceLayout::Pooled>>;
using UnpooledSpikeSlabPriorAD = IndependentTopology<
    mode_ad,
    SpikeSlabPrior<VarianceLayout::Unpooled>,
    SpikeSlabPrior<VarianceLayout::Unpooled>>;
using FixedUnpooledSpikeSlabPriorAD = IndependentTopology<
    mode_ad,
    SpikeSlabPrior<VarianceLayout::Unpooled, UpdatePolicy::Fixed>,
    SpikeSlabPrior<VarianceLayout::Unpooled, UpdatePolicy::Fixed>>;
using HeterogeneousPriorAD = IndependentTopology<
    mode_ad,
    GaussianPrior<VarianceLayout::Pooled>,
    SpikeSlabPrior<VarianceLayout::Unpooled, UpdatePolicy::Fixed>>;
using ScaledMixturePriorAD = IndependentTopology<
    mode_ad,
    ScaledMixturePrior<ScaledMixture::class_count>,
    ScaledMixturePrior<ScaledMixture::class_count>>;
using JointPrior = JointTopology<
    IndependentTopology<
        mode_ad,
        GaussianPrior<VarianceLayout::Pooled>,
        HalfNormalPrior<HalfNormalAsymmetry::Count>>,
    JointSpikeSlabPrior<JointSpikeSlab::class_count>>;
using MagnitudeJointPrior = JointTopology<
    IndependentTopology<
        mode_ad,
        GaussianPrior<VarianceLayout::Pooled>,
        HalfNormalPrior<HalfNormalAsymmetry::Magnitude>>,
    JointSpikeSlabPrior<JointSpikeSlab::class_count>>;
using JointModeParameters = IndependentTopology<
    mode_ad,
    Gaussian,
    HalfNormal<HalfNormalAsymmetry::Count>>;
using JointSpikeSlabAD = JointTopology<JointModeParameters, JointSpikeSlab>;

template <std::size_t ClassCount>
concept ValidScaledMixtureState
    = requires { typename ScaledMixtureState<ClassCount>; };

template <std::size_t ClassCount>
concept ValidJointSpikeSlabState
    = requires { typename JointSpikeSlabState<ClassCount>; };

static_assert(
    ScaledMixtureState<ScaledMixture::class_count>::component_count == 4);
static_assert(
    JointSpikeSlabState<JointSpikeSlab::class_count>::component_count == 3);
static_assert(ValidScaledMixtureState<256>);
static_assert(!ValidScaledMixtureState<257>);
static_assert(ValidJointSpikeSlabState<256>);
static_assert(!ValidJointSpikeSlabState<257>);

static_assert(std::same_as<
              decltype(SpikeSlabState<VarianceLayout::Pooled>::assignment),
              Eigen::VectorX<std::uint8_t>>);
static_assert(
    std::same_as<
        decltype(ScaledMixtureState<ScaledMixture::class_count>::assignment),
        Eigen::VectorX<std::uint8_t>>);
static_assert(
    std::same_as<
        decltype(JointSpikeSlabState<JointSpikeSlab::class_count>::assignment),
        Eigen::VectorX<std::uint8_t>>);
static_assert(std::same_as<
              decltype(HalfNormalState<HalfNormalAsymmetry::Count>::assignment),
              Eigen::VectorX<std::uint8_t>>);

static_assert(std::same_as<
              gelex::detail::genetic_state_t<PooledGaussianPriorAD>,
              IndependentTopology<
                  mode_ad,
                  GeneticModeState<GaussianState<VarianceLayout::Pooled>>,
                  GeneticModeState<GaussianState<VarianceLayout::Pooled>>>>);
static_assert(std::same_as<
              gelex::detail::genetic_state_t<UnpooledGaussianPriorA>,
              IndependentTopology<
                  mode_a,
                  GeneticModeState<GaussianState<VarianceLayout::Unpooled>>>>);
static_assert(std::same_as<
              gelex::detail::genetic_state_t<PooledSpikeSlabPriorAD>,
              IndependentTopology<
                  mode_ad,
                  GeneticModeState<SpikeSlabState<VarianceLayout::Pooled>>,
                  GeneticModeState<SpikeSlabState<VarianceLayout::Pooled>>>>);
static_assert(std::same_as<
              gelex::detail::genetic_state_t<UnpooledSpikeSlabPriorAD>,
              IndependentTopology<
                  mode_ad,
                  GeneticModeState<SpikeSlabState<VarianceLayout::Unpooled>>,
                  GeneticModeState<SpikeSlabState<VarianceLayout::Unpooled>>>>);
static_assert(std::same_as<
              gelex::detail::genetic_state_t<UnpooledSpikeSlabPriorAD>,
              gelex::detail::genetic_state_t<FixedUnpooledSpikeSlabPriorAD>>);
static_assert(std::same_as<
              gelex::detail::genetic_state_t<HeterogeneousPriorAD>,
              IndependentTopology<
                  mode_ad,
                  GeneticModeState<GaussianState<VarianceLayout::Pooled>>,
                  GeneticModeState<SpikeSlabState<VarianceLayout::Unpooled>>>>);
static_assert(
    std::same_as<
        gelex::detail::genetic_state_t<ScaledMixturePriorAD>,
        IndependentTopology<
            mode_ad,
            GeneticModeState<ScaledMixtureState<ScaledMixture::class_count>>,
            GeneticModeState<ScaledMixtureState<ScaledMixture::class_count>>>>);
static_assert(
    std::same_as<
        gelex::detail::genetic_state_t<JointPrior>,
        JointTopology<
            IndependentTopology<
                mode_ad,
                GeneticModeState<GaussianState<VarianceLayout::Pooled>>,
                GeneticModeState<HalfNormalState<HalfNormalAsymmetry::Count>>>,
            JointSpikeSlabState<JointSpikeSlab::class_count>>>);
static_assert(std::same_as<
              gelex::detail::genetic_state_t<MagnitudeJointPrior>,
              JointTopology<
                  IndependentTopology<
                      mode_ad,
                      GeneticModeState<GaussianState<VarianceLayout::Pooled>>,
                      GeneticModeState<
                          HalfNormalState<HalfNormalAsymmetry::Magnitude>>>,
                  JointSpikeSlabState<JointSpikeSlab::class_count>>>);
static_assert(std::same_as<
              gelex::bayes_state_t<BayesPrior<PooledGaussianPriorA>>,
              BayesState<PooledGaussianPriorA>>);

using PooledGaussianMethod = GaussianMethod<VarianceLayout::Pooled>;
using UnpooledGaussianMethod = GaussianMethod<VarianceLayout::Unpooled>;
using FixedUnpooledSpikeSlabMethod
    = SpikeSlabMethod<VarianceLayout::Unpooled, UpdatePolicy::Fixed>;
using DefaultScaledMixtureMethod = ScaledMixtureMethod<>;
using FixedJointSpikeSlabMethod = JointSpikeSlabMethod<UpdatePolicy::Fixed>;
using MagnitudeJointSpikeSlabMethod = JointSpikeSlabMethod<
    UpdatePolicy::Sampled,
    HalfNormalAsymmetry::Magnitude>;
using SpikeSlabAD = IndependentTopology<mode_ad, SpikeSlab, SpikeSlab>;

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

    state.for_each(
        [&]<GeneticMode Mode>(const auto& mode_state)
        {
            STATIC_REQUIRE(mode_ad.contains(Mode));
            REQUIRE(mode_state.coefficients.size() == model.genetic().cols());
            REQUIRE(mode_state.coefficients.isZero());
            REQUIRE(
                mode_state.family_state.fitted_values.size()
                == model.genetic().rows());
            REQUIRE(mode_state.family_state.fitted_values.isZero());
        });
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
            SpikeSlab{0.05},
            SpikeSlab{0.2},
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
    REQUIRE(family_state.fitted_values.rows() == model.genetic().rows());
    REQUIRE(
        family_state.fitted_values.cols()
        == static_cast<Eigen::Index>(
            ScaledMixtureState<ScaledMixture::class_count>::component_count));
    REQUIRE(family_state.fitted_values.isZero());
}

TEST_CASE(
    "joint state keeps shared latent variables in the joint value",
    "[bayes][state_compilation]")
{
    const auto model = make_model(mode_ad);
    const auto recipe = BayesRecipe<FixedJointSpikeSlabMethod, mode_ad>{
        JointSpikeSlabAD{
            JointModeParameters{
                Gaussian{}, HalfNormal<HalfNormalAsymmetry::Count>{0.6}},
            JointSpikeSlab{{0.8, 0.1, 0.05, 0.05}}},
        VarianceBudget{{.additive = 0.4, .dominance = 0.1}}};
    const auto prior = gelex::compile(recipe, model);

    const auto state
        = gelex::detail::make_genetic_state(prior.genetic(), model.genetic());
    const auto& joint = state.joint();

    REQUIRE(joint.probabilities == recipe.parameters().joint().probabilities());
    REQUIRE(joint.assignment.size() == model.genetic().cols());
    REQUIRE(joint.assignment.isZero());
    REQUIRE(joint.fitted_values.rows() == model.genetic().rows());
    REQUIRE(
        joint.fitted_values.cols()
        == static_cast<Eigen::Index>(
            JointSpikeSlabState<JointSpikeSlab::class_count>::component_count));
    REQUIRE(joint.fitted_values.isZero());
    state.mode_values().for_each(
        [&]<GeneticMode Mode>(const auto& mode_state)
        {
            STATIC_REQUIRE(mode_ad.contains(Mode));
            REQUIRE(
                mode_state.family_state.fitted_values.size()
                == model.genetic().rows());
            REQUIRE(mode_state.family_state.fitted_values.isZero());
        });
    const auto& dominance
        = state.mode_values().get<GeneticMode::D>().family_state;
    REQUIRE(dominance.positive_probability == 0.6);
    REQUIRE(dominance.assignment.size() == model.genetic().cols());
    REQUIRE(dominance.assignment.isZero());
    REQUIRE(
        state.mode_values().get<GeneticMode::A>().family_state.variance
        == Approx(prior.genetic()
                      .mode_values()
                      .get<GeneticMode::A>()
                      .variance.initial_value()));
}

TEST_CASE(
    "magnitude half-normal state initializes both variance axes",
    "[bayes][state_compilation]")
{
    const auto model = make_model(mode_ad);
    const auto prior = gelex::compile(
        BayesRecipe<MagnitudeJointSpikeSlabMethod, mode_ad>::defaults(), model);

    const auto state
        = gelex::detail::make_genetic_state(prior.genetic(), model.genetic());
    const auto& dominance_prior
        = prior.genetic().mode_values().get<GeneticMode::D>();
    const auto& dominance_state
        = state.mode_values().get<GeneticMode::D>().family_state;
    const double initial_variance = dominance_prior.variance.initial_value();

    REQUIRE(
        dominance_state.variances
        == std::array{initial_variance, initial_variance});
    REQUIRE(dominance_state.assignment.size() == model.genetic().cols());
    REQUIRE(dominance_state.assignment.isZero());
    REQUIRE(dominance_state.fitted_values.size() == model.genetic().rows());
    REQUIRE(dominance_state.fitted_values.isZero());
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
    REQUIRE(genetic.family_state.fitted_values.isZero());
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
