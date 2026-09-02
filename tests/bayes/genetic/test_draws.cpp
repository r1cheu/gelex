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

#include "gelex/bayes/detail/draws_factory.h"
#include "gelex/bayes/detail/state_factory.h"
#include "gelex/bayes/draw.h"
#include "gelex/bayes/genetic/draws.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/spec.h"
#include "gelex/io/binary_reader.h"
#include "gelex/io/binary_writer.h"
#include "gelex/types/genetic_mode.h"
#include "gelex/types/mode_values.h"

#include "compact_genotype_fixture.h"
#include "file_fixture.h"

using Catch::Approx;
using gelex::BayesModel;
using gelex::BayesRecipe;
using gelex::EmptyDraw;
using gelex::GaussianDraws;
using gelex::GaussianFamily;
using gelex::GaussianPrior;
using gelex::GeneticMode;
using gelex::GeneticModeDraws;
using gelex::GeneticModeSet;
using gelex::HalfNormalAsymmetry;
using gelex::HalfNormalDraws;
using gelex::HalfNormalPrior;
using gelex::IndependentDraws;
using gelex::JointDraws;
using gelex::JointModeValues;
using gelex::JointSpikeSlab;
using gelex::JointSpikeSlabDraws;
using gelex::JointSpikeSlabFamily;
using gelex::JointSpikeSlabPrior;
using gelex::ModeValues;
using gelex::ScalarDraw;
using gelex::ScaledMixture;
using gelex::ScaledMixtureDraws;
using gelex::ScaledMixtureFamily;
using gelex::SpikeSlabDraws;
using gelex::SpikeSlabFamily;
using gelex::SpikeSlabPrior;
using gelex::UpdatePolicy;
using gelex::VarianceLayout;
using gelex::detail::genetic_draws_t;
using gelex::detail::genetic_state_t;

namespace
{

constexpr auto mode_a = GeneticModeSet{GeneticMode::A};
constexpr auto mode_ad = GeneticMode::A | GeneticMode::D;

using PooledGaussianPriorAD = ModeValues<
    mode_ad,
    GaussianPrior<VarianceLayout::Pooled>,
    GaussianPrior<VarianceLayout::Pooled>>;
using UnpooledSpikeSlabPriorA
    = ModeValues<mode_a, SpikeSlabPrior<VarianceLayout::Unpooled>>;
using FixedUnpooledSpikeSlabPriorA = ModeValues<
    mode_a,
    SpikeSlabPrior<VarianceLayout::Unpooled, UpdatePolicy::Fixed>>;
using JointPrior = JointModeValues<
    ModeValues<
        mode_ad,
        GaussianPrior<VarianceLayout::Pooled>,
        HalfNormalPrior<HalfNormalAsymmetry::Count>>,
    JointSpikeSlabPrior<JointSpikeSlab::class_count>>;

using UnpooledSpikeSlabFamily = SpikeSlabFamily<VarianceLayout::Unpooled>;
using FixedPooledSpikeSlabFamily
    = SpikeSlabFamily<VarianceLayout::Pooled, UpdatePolicy::Fixed>;
using MagnitudeJointSpikeSlabFamily = JointSpikeSlabFamily<
    UpdatePolicy::Sampled,
    HalfNormalAsymmetry::Magnitude>;

static_assert(
    std::same_as<
        genetic_draws_t<PooledGaussianPriorAD>,
        IndependentDraws<
            genetic_state_t<PooledGaussianPriorAD>,
            ModeValues<
                mode_ad,
                GeneticModeDraws<GaussianDraws<VarianceLayout::Pooled>>,
                GeneticModeDraws<GaussianDraws<VarianceLayout::Pooled>>>>>);
static_assert(std::same_as<
              genetic_draws_t<UnpooledSpikeSlabPriorA>,
              IndependentDraws<
                  genetic_state_t<UnpooledSpikeSlabPriorA>,
                  ModeValues<
                      mode_a,
                      GeneticModeDraws<SpikeSlabDraws<
                          VarianceLayout::Unpooled,
                          UpdatePolicy::Sampled>>>>>);
static_assert(
    std::same_as<
        genetic_draws_t<JointPrior>,
        JointDraws<
            genetic_state_t<JointPrior>,
            ModeValues<
                mode_ad,
                GeneticModeDraws<GaussianDraws<VarianceLayout::Pooled>>,
                GeneticModeDraws<HalfNormalDraws<HalfNormalAsymmetry::Count>>>,
            JointSpikeSlabDraws<
                JointSpikeSlab::class_count,
                UpdatePolicy::Sampled>>>);

// Unlike state, the draws tree does distinguish fixed from sampled parameters:
// a fixed parameter maps to EmptyDraw and reserves no payload.
static_assert(std::same_as<
              genetic_state_t<UnpooledSpikeSlabPriorA>,
              genetic_state_t<FixedUnpooledSpikeSlabPriorA>>);
static_assert(!std::same_as<
              genetic_draws_t<UnpooledSpikeSlabPriorA>,
              genetic_draws_t<FixedUnpooledSpikeSlabPriorA>>);
static_assert(std::same_as<
              decltype(SpikeSlabDraws<
                       VarianceLayout::Unpooled,
                       UpdatePolicy::Sampled>::probability),
              ScalarDraw>);
static_assert(
    std::same_as<
        decltype(SpikeSlabDraws<VarianceLayout::Unpooled, UpdatePolicy::Fixed>::
                     probability),
        EmptyDraw>);
static_assert(std::same_as<
              decltype(ScaledMixtureDraws<
                       ScaledMixture::class_count,
                       UpdatePolicy::Fixed>::probabilities),
              EmptyDraw>);

auto make_model(GeneticModeSet modes) -> BayesModel
{
    return gelex::test::make_compact_model(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}},
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        modes);
}

}  // namespace

TEST_CASE(
    "unpooled spike-slab draws record every sampled leaf parameter",
    "[bayes][draws][genetic]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "spike_slab.draws";
    const auto model = make_model(mode_a);
    const auto prior = gelex::make_prior(
        BayesRecipe<mode_a, UnpooledSpikeSlabFamily>::defaults(), model);
    auto state
        = gelex::detail::make_genetic_state(prior.genetic(), model.genetic());

    {
        gelex::BinaryWriter writer(path.string());
        auto draws = gelex::detail::make_genetic_draws(
            prior.genetic(), model.genetic(), writer, 2);

        auto& additive = state.get<GeneticMode::A>();
        additive.coefficients = Eigen::VectorXd{{1.0, 2.0}};
        additive.family_state.variance = Eigen::VectorXd{{0.1, 0.2}};
        additive.family_state.assignment = Eigen::VectorX<std::uint8_t>{{0, 1}};
        additive.family_state.probability = 0.3;
        draws.append(state);

        additive.coefficients = Eigen::VectorXd{{3.0, 4.0}};
        additive.family_state.variance = Eigen::VectorXd{{0.3, 0.4}};
        additive.family_state.assignment = Eigen::VectorX<std::uint8_t>{{1, 1}};
        additive.family_state.probability = 0.5;
        draws.append(state);
        writer.close();

        const auto& leaf = draws.get<GeneticMode::A>();
        REQUIRE(leaf.coefficients.result().mean.isApprox(
            Eigen::VectorXd{{2.0, 3.0}}));
        REQUIRE(leaf.family_draws.probability.result().mean == Approx(0.4));
        REQUIRE(leaf.family_draws.assignment.result().probabilities.isApprox(
            Eigen::MatrixXd{{0.5, 0.5}, {0.0, 1.0}}));
    }

    const gelex::BinaryReader reader(path.string());
    REQUIRE(reader.to_map<float>("genetic/A/coefficients")
                .isApprox(Eigen::MatrixXf{{1.0, 3.0}, {2.0, 4.0}}));
    REQUIRE(reader.to_map<float>("genetic/A/variance")
                .isApprox(Eigen::MatrixXf{{0.1, 0.3}, {0.2, 0.4}}));
    REQUIRE(reader.to_map<double>("genetic/A/probability")
                .isApprox(Eigen::MatrixXd{{0.3, 0.5}}));
    REQUIRE(reader.to_map<std::uint8_t>("genetic/A/assignment")
                .isApprox(
                    Eigen::Matrix<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic>{
                        {0, 1}, {1, 1}}));
}

TEST_CASE("a fixed probability reserves no payload", "[bayes][draws][genetic]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "fixed_spike_slab.draws";
    const auto model = make_model(mode_a);
    const auto prior = gelex::make_prior(
        BayesRecipe<mode_a, FixedPooledSpikeSlabFamily>::defaults(), model);
    auto state
        = gelex::detail::make_genetic_state(prior.genetic(), model.genetic());

    {
        gelex::BinaryWriter writer(path.string());
        auto draws = gelex::detail::make_genetic_draws(
            prior.genetic(), model.genetic(), writer, 1);
        draws.append(state);
        writer.close();
    }

    const gelex::BinaryReader reader(path.string());
    REQUIRE(reader.contains("genetic/A/coefficients"));
    REQUIRE(reader.contains("genetic/A/variance"));
    REQUIRE(reader.contains("genetic/A/assignment"));
    REQUIRE_FALSE(reader.contains("genetic/A/probability"));
}

TEST_CASE(
    "scaled-mixture draws record the class simplex per draw",
    "[bayes][draws][genetic]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "scaled_mixture.draws";
    const auto model = make_model(mode_a);
    const auto prior = gelex::make_prior(
        BayesRecipe<mode_a, ScaledMixtureFamily<>>::defaults(), model);
    auto state
        = gelex::detail::make_genetic_state(prior.genetic(), model.genetic());

    {
        gelex::BinaryWriter writer(path.string());
        auto draws = gelex::detail::make_genetic_draws(
            prior.genetic(), model.genetic(), writer, 1);

        auto& additive = state.get<GeneticMode::A>();
        additive.family_state.assignment = Eigen::VectorX<std::uint8_t>{{0, 4}};
        additive.family_state.probabilities = {0.5, 0.2, 0.15, 0.1, 0.05};
        draws.append(state);
        writer.close();
    }

    const gelex::BinaryReader reader(path.string());
    REQUIRE(
        reader.to_map<float>("genetic/A/probabilities")
            .isApprox(Eigen::MatrixXf{{0.5}, {0.2}, {0.15}, {0.1}, {0.05}}));
    REQUIRE(reader.to_map<std::uint8_t>("genetic/A/assignment")
                .isApprox(
                    Eigen::Matrix<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic>{
                        {0}, {4}}));
}

TEST_CASE(
    "joint spike-slab draws record both mode leaves and the joint leaf",
    "[bayes][draws][genetic]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "joint.draws";
    const auto model = make_model(mode_ad);
    const auto prior = gelex::make_prior(
        BayesRecipe<mode_ad, JointSpikeSlabFamily<>>::defaults(), model);
    auto state
        = gelex::detail::make_genetic_state(prior.genetic(), model.genetic());

    {
        gelex::BinaryWriter writer(path.string());
        auto draws = gelex::detail::make_genetic_draws(
            prior.genetic(), model.genetic(), writer, 1);

        auto& additive = state.mode_values().get<GeneticMode::A>();
        additive.coefficients = Eigen::VectorXd{{1.0, 2.0}};
        additive.family_state.variance = 0.5;

        auto& dominance = state.mode_values().get<GeneticMode::D>();
        dominance.coefficients = Eigen::VectorXd{{3.0, 4.0}};
        dominance.family_state.variance = 0.25;
        dominance.family_state.assignment
            = Eigen::VectorX<std::uint8_t>{{0, 2}};
        dominance.family_state.positive_probability = 0.7;

        state.joint().assignment = Eigen::VectorX<std::uint8_t>{{1, 3}};
        state.joint().probabilities = {0.7, 0.1, 0.1, 0.1};
        draws.append(state);
        writer.close();

        REQUIRE(
            draws.get<GeneticMode::D>()
                .family_draws.positive_probability.result()
                .mean
            == Approx(0.7));
        REQUIRE(draws.get<GeneticMode::D>()
                    .family_draws.assignment.result()
                    .probabilities.isApprox(
                        Eigen::MatrixXd{{1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}}));
    }

    const gelex::BinaryReader reader(path.string());
    REQUIRE(reader.to_map<float>("genetic/A/coefficients")
                .isApprox(Eigen::MatrixXf{{1.0}, {2.0}}));
    REQUIRE(reader.to_map<double>("genetic/A/variance")
                .isApprox(Eigen::MatrixXd{{0.5}}));
    REQUIRE(reader.to_map<double>("genetic/D/variance")
                .isApprox(Eigen::MatrixXd{{0.25}}));
    REQUIRE(reader.to_map<double>("genetic/D/positive_probability")
                .isApprox(Eigen::MatrixXd{{0.7}}));
    REQUIRE(reader.to_map<std::uint8_t>("genetic/D/assignment")
                .isApprox(
                    Eigen::Matrix<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic>{
                        {0}, {2}}));
    REQUIRE(reader.to_map<std::uint8_t>("genetic/joint/assignment")
                .isApprox(
                    Eigen::Matrix<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic>{
                        {1}, {3}}));
    REQUIRE(reader.to_map<float>("genetic/joint/probabilities")
                .isApprox(Eigen::MatrixXf{{0.7}, {0.1}, {0.1}, {0.1}}));
}

TEST_CASE(
    "magnitude half-normal draws record both variance axes",
    "[bayes][draws][genetic]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "magnitude.draws";
    const auto model = make_model(mode_ad);
    const auto prior = gelex::make_prior(
        BayesRecipe<mode_ad, MagnitudeJointSpikeSlabFamily>::defaults(), model);
    auto state
        = gelex::detail::make_genetic_state(prior.genetic(), model.genetic());

    {
        gelex::BinaryWriter writer(path.string());
        auto draws = gelex::detail::make_genetic_draws(
            prior.genetic(), model.genetic(), writer, 1);

        auto& dominance = state.mode_values().get<GeneticMode::D>();
        dominance.family_state.variances = {0.2, 0.8};
        dominance.family_state.assignment
            = Eigen::VectorX<std::uint8_t>{{2, 0}};
        draws.append(state);
        writer.close();

        REQUIRE(draws.get<GeneticMode::D>()
                    .family_draws.assignment.result()
                    .probabilities.isApprox(
                        Eigen::MatrixXd{{0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}}));
    }

    const gelex::BinaryReader reader(path.string());
    REQUIRE(reader.to_map<float>("genetic/D/variance")
                .isApprox(Eigen::MatrixXf{{0.2}, {0.8}}));
    REQUIRE(reader.to_map<std::uint8_t>("genetic/D/assignment")
                .isApprox(
                    Eigen::Matrix<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic>{
                        {2}, {0}}));
    REQUIRE_FALSE(reader.contains("genetic/D/positive_probability"));
}
