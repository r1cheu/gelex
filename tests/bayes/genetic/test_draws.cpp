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

#include "gelex/bayes/basic_draw.h"
#include "gelex/bayes/detail/draws_factory.h"
#include "gelex/bayes/detail/pip_factory.h"
#include "gelex/bayes/detail/state_factory.h"
#include "gelex/bayes/genetic/draws.h"
#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/spike_slab.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/spec.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_reader.h"
#include "gelex/io/binary_writer.h"

#include "compact_genotype_fixture.h"
#include "file_fixture.h"

using Catch::Approx;
using gelex::BayesModel;
using gelex::BayesRecipe;
using gelex::EmptyDraw;
using gelex::GaussianDraws;
using gelex::GaussianFamily;
using gelex::GaussianPrior;
using gelex::GeneticCoefficientDraws;
using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::HalfNormalDraws;
using gelex::HalfNormalPrior;
using gelex::IndependentGeneticDraws;
using gelex::JointGeneticDraws;
using gelex::JointModeValues;
using gelex::JointSpikeSlab;
using gelex::JointSpikeSlabDraws;
using gelex::JointSpikeSlabFamily;
using gelex::JointSpikeSlabPrior;
using gelex::MixtureWeightUpdate;
using gelex::ModeValues;
using gelex::ScalarDraw;
using gelex::ScaledMixture;
using gelex::ScaledMixtureDraws;
using gelex::ScaledMixtureFamily;
using gelex::SpikeSlabDraws;
using gelex::SpikeSlabFamily;
using gelex::SpikeSlabPrior;
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
    SpikeSlabPrior<VarianceLayout::Unpooled, MixtureWeightUpdate::Disabled>>;
using JointPrior = JointModeValues<
    ModeValues<mode_ad, GaussianPrior<VarianceLayout::Pooled>, HalfNormalPrior>,
    JointSpikeSlabPrior<JointSpikeSlab::class_count>>;

using UnpooledSpikeSlabFamily = SpikeSlabFamily<VarianceLayout::Unpooled>;
using FixedPooledSpikeSlabFamily
    = SpikeSlabFamily<VarianceLayout::Pooled, MixtureWeightUpdate::Disabled>;

static_assert(std::same_as<
              genetic_draws_t<PooledGaussianPriorAD>,
              IndependentGeneticDraws<
                  GeneticCoefficientDraws<mode_ad>,
                  ModeValues<
                      mode_ad,
                      GaussianDraws<VarianceLayout::Pooled>,
                      GaussianDraws<VarianceLayout::Pooled>>>>);
static_assert(std::same_as<
              genetic_draws_t<UnpooledSpikeSlabPriorA>,
              IndependentGeneticDraws<
                  GeneticCoefficientDraws<mode_a>,
                  ModeValues<
                      mode_a,
                      SpikeSlabDraws<
                          VarianceLayout::Unpooled,
                          MixtureWeightUpdate::Enabled>>>>);
static_assert(std::same_as<
              genetic_draws_t<JointPrior>,
              JointGeneticDraws<
                  GeneticCoefficientDraws<mode_ad>,
                  ModeValues<
                      mode_ad,
                      GaussianDraws<VarianceLayout::Pooled>,
                      HalfNormalDraws>,
                  JointSpikeSlabDraws<
                      JointSpikeSlab::class_count,
                      MixtureWeightUpdate::Enabled>>>);

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
                       MixtureWeightUpdate::Enabled>::probability),
              ScalarDraw>);
static_assert(std::same_as<
              decltype(SpikeSlabDraws<
                       VarianceLayout::Unpooled,
                       MixtureWeightUpdate::Disabled>::probability),
              EmptyDraw>);
static_assert(std::same_as<
              decltype(ScaledMixtureDraws<
                       ScaledMixture::class_count,
                       MixtureWeightUpdate::Disabled>::probabilities),
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
    auto state = gelex::detail::make_state(prior.genetic(), model.genetic());

    {
        gelex::BinaryWriter writer(path.string());
        auto draws = gelex::detail::make_draws(
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

        REQUIRE(
            draws.coefficients().get<GeneticMode::A>().result().mean.isApprox(
                Eigen::VectorXd{{2.0, 3.0}}));
        const auto& family = draws.family<GeneticMode::A>();
        REQUIRE(family.probability.result().mean == Approx(0.4));
        REQUIRE(family.assignment.result().probabilities.isApprox(
            Eigen::MatrixXd{{0.5, 0.5}, {0.0, 1.0}}));
        REQUIRE(
            gelex::detail::make_pip(draws)
                .get<GeneticMode::A>()
                .probabilities()
                .isApprox(Eigen::VectorXd{{0.5, 1.0}}));
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
    auto state = gelex::detail::make_state(prior.genetic(), model.genetic());

    {
        gelex::BinaryWriter writer(path.string());
        auto draws = gelex::detail::make_draws(
            prior.genetic(), model.genetic(), writer, 1);
        draws.append(state);
    }

    const gelex::BinaryReader reader(path.string());
    REQUIRE(reader.contains("genetic/A/coefficients"));
    REQUIRE(reader.contains("genetic/A/variance"));
    REQUIRE(reader.contains("genetic/A/assignment"));
    REQUIRE_FALSE(reader.contains("genetic/A/probability"));
}

TEST_CASE(
    "scaled-mixture draws record the class simplex and per-class explained "
    "variance",
    "[bayes][draws][genetic]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "scaled_mixture.draws";
    const auto model = make_model(mode_a);
    const auto prior = gelex::make_prior(
        BayesRecipe<mode_a, ScaledMixtureFamily<>>::defaults(), model);
    auto state = gelex::detail::make_state(prior.genetic(), model.genetic());

    {
        gelex::BinaryWriter writer(path.string());
        auto draws = gelex::detail::make_draws(
            prior.genetic(), model.genetic(), writer, 1);

        auto& additive = state.get<GeneticMode::A>();
        additive.family_state.assignment = Eigen::VectorX<std::uint8_t>{{0, 4}};
        additive.family_state.probabilities = {0.5, 0.2, 0.15, 0.1, 0.05};
        additive.family_state.fitted_values = Eigen::MatrixXd{
            {0.0, 1.0, 2.0, 0.0}, {3.0, 1.0, 0.0, 0.0}, {0.0, 1.0, 4.0, 0.0}};
        draws.append(state);

        REQUIRE(
            gelex::detail::make_pip(draws)
                .get<GeneticMode::A>()
                .probabilities()
                .isApprox(Eigen::VectorXd{{0.0, 1.0}}));
    }

    const gelex::BinaryReader reader(path.string());
    REQUIRE(
        reader.to_map<float>("genetic/A/probabilities")
            .isApprox(Eigen::MatrixXf{{0.5}, {0.2}, {0.15}, {0.1}, {0.05}}));
    REQUIRE(reader.to_map<std::uint8_t>("genetic/A/assignment")
                .isApprox(
                    Eigen::Matrix<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic>{
                        {0}, {4}}));
    REQUIRE(reader.to_map<float>("genetic/A/component_explained_variance")
                .isApprox(Eigen::MatrixXf{{2.0}, {0.0}, {8.0 / 3.0}, {0.0}}));
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
    auto state = gelex::detail::make_state(prior.genetic(), model.genetic());

    {
        gelex::BinaryWriter writer(path.string());
        auto draws = gelex::detail::make_draws(
            prior.genetic(), model.genetic(), writer, 1);

        auto& additive = state.mode_values().get<GeneticMode::A>();
        additive.coefficients = Eigen::VectorXd{{1.0, 2.0}};
        additive.family_state.variance = 0.5;

        auto& dominance = state.mode_values().get<GeneticMode::D>();
        dominance.coefficients = Eigen::VectorXd{{3.0, 4.0}};
        dominance.family_state.variance = 0.25;
        dominance.family_state.probit_coefficients
            = Eigen::Vector2d{{0.7, -0.2}};

        state.joint().assignment = Eigen::VectorX<std::uint8_t>{{1, 3}};
        state.joint().probabilities = {0.7, 0.1, 0.1, 0.1};
        // Columns are A in A-only, A in AD, D in D-only, D in AD.
        state.joint().fitted_values = Eigen::MatrixXd{
            {1.0, 0.0, 2.0, 0.0}, {2.0, 0.0, 2.0, 4.0}, {3.0, 3.0, 2.0, 2.0}};
        draws.append(state);

        REQUIRE(draws.family<GeneticMode::D>()
                    .probit_coefficients.result()
                    .mean.isApprox(Eigen::Vector2d{{0.7, -0.2}}));
        const auto pip = gelex::detail::make_pip(draws);
        REQUIRE(pip.get<GeneticMode::A>().probabilities().isApprox(
            Eigen::VectorXd{{1.0, 1.0}}));
        REQUIRE(pip.get<GeneticMode::D>().probabilities().isApprox(
            Eigen::VectorXd{{0.0, 1.0}}));
        REQUIRE(
            pip.joint().probabilities().isApprox(Eigen::VectorXd{{1.0, 1.0}}));
    }

    const gelex::BinaryReader reader(path.string());
    REQUIRE(reader.to_map<float>("genetic/A/coefficients")
                .isApprox(Eigen::MatrixXf{{1.0}, {2.0}}));
    REQUIRE(reader.to_map<double>("genetic/A/variance")
                .isApprox(Eigen::MatrixXd{{0.5}}));
    REQUIRE(reader.to_map<double>("genetic/D/variance")
                .isApprox(Eigen::MatrixXd{{0.25}}));
    REQUIRE(reader.to_map<float>("genetic/D/probit_coefficients")
                .isApprox(Eigen::MatrixXf{{0.7F}, {-0.2F}}));
    REQUIRE(reader.to_map<std::uint8_t>("genetic/joint/assignment")
                .isApprox(
                    Eigen::Matrix<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic>{
                        {1}, {3}}));
    REQUIRE(reader.to_map<float>("genetic/joint/probabilities")
                .isApprox(Eigen::MatrixXf{{0.7}, {0.1}, {0.1}, {0.1}}));
    REQUIRE(
        reader.to_map<float>("genetic/joint/component_explained_variance")
            .isApprox(Eigen::MatrixXf{{2.0 / 3.0}, {2.0}, {0.0}, {8.0 / 3.0}}));
}
