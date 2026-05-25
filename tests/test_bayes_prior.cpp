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

#include <array>
#include <memory>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/exception.h"
#include "gelex/model/bayes/capabilities.h"
#include "gelex/model/bayes/gaussian_prior.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/prior_state.h"
#include "gelex/model/bayes/state_capabilities.h"
#include "gelex/types/genetic_effect_type.h"

using gelex::GelexException;
using gelex::GeneticMode;
using gelex::bayes::BayesPrior;
using gelex::bayes::ComponentStateCap;
using gelex::bayes::DirichletPrior;
using gelex::bayes::GaussianPrior;
using gelex::bayes::GeneticPrior;
using gelex::bayes::JointMixtureGaussianPrior;
using gelex::bayes::MarkerVariance;
using gelex::bayes::MarkerVarianceCap;
using gelex::bayes::MarkerVarianceLayout;
using gelex::bayes::MixtureProportion;
using gelex::bayes::MixtureProportionCap;
using gelex::bayes::MultiplierCap;
using gelex::bayes::ProportionStateCap;
using gelex::bayes::RandomPrior;
using gelex::bayes::ResidualPrior;
using gelex::bayes::ScaledInvChiSqPrior;
using gelex::bayes::ScaledMixtureGaussianPrior;
using gelex::bayes::SimplexParameter;
using gelex::bayes::SpikeSlabGaussianPrior;
using gelex::bayes::UpdatePolicy;
using gelex::bayes::VarianceParameter;
using gelex::bayes::VarianceStateCap;

namespace
{

auto make_variance(double init = 1.0) -> VarianceParameter
{
    return VarianceParameter(init, ScaledInvChiSqPrior{4.0, 1.0});
}

auto make_random_prior(double init = 1.0) -> RandomPrior
{
    return RandomPrior{make_variance(init)};
}

auto make_residual_prior(double init = 1.0) -> ResidualPrior
{
    return ResidualPrior{make_variance(init)};
}

auto make_marker_variance() -> MarkerVariance
{
    return MarkerVariance{
        MarkerVarianceLayout::per_marker,
        make_variance(),
    };
}

auto make_proportion_2() -> MixtureProportion
{
    return MixtureProportion{
        SimplexParameter{
            Eigen::VectorXd{{0.9, 0.1}},
            DirichletPrior{Eigen::VectorXd{{1.0, 1.0}}}},
        UpdatePolicy::fixed,
    };
}

auto make_multiplier_2() -> Eigen::VectorXd
{
    return Eigen::VectorXd{{0.0, 1.0}};
}

auto make_gaussian(GeneticMode mode = GeneticMode::A)
    -> std::unique_ptr<GeneticPrior>
{
    return std::make_unique<GaussianPrior>(mode, make_marker_variance());
}

auto make_spike_slab(GeneticMode mode = GeneticMode::A)
    -> std::unique_ptr<GeneticPrior>
{
    return std::make_unique<SpikeSlabGaussianPrior>(
        mode, make_marker_variance(), make_proportion_2());
}

auto make_joint() -> std::unique_ptr<GeneticPrior>
{
    return std::make_unique<JointMixtureGaussianPrior>(
        std::array<GeneticMode, 2>{GeneticMode::A, GeneticMode::D},
        std::array<MarkerVariance, 2>{
            make_marker_variance(), make_marker_variance()},
        make_proportion_2());
}

}  // namespace

TEST_CASE("BayesPrior accessors expose construction arguments", "[bayes_prior]")
{
    std::vector<std::unique_ptr<GeneticPrior>> genetics;
    genetics.push_back(make_gaussian(GeneticMode::A));
    genetics.push_back(make_spike_slab(GeneticMode::D));

    BayesPrior prior(
        make_random_prior(2.0), std::move(genetics), make_residual_prior(5.0));

    REQUIRE(prior.random().initial_value() == 2.0);
    REQUIRE(prior.residual().initial_value() == 5.0);

    const auto range = prior.genetics();
    REQUIRE(range.size() == 2);
    REQUIRE_FALSE(range.empty());

    auto it = range.begin();
    REQUIRE((*it).contains(GeneticMode::A));
    ++it;
    REQUIRE((*it).contains(GeneticMode::D));
    ++it;
    REQUIRE(it == range.end());
}

TEST_CASE("GeneticPrior capabilities compose prior data axes", "[bayes_prior]")
{
    GaussianPrior gaussian(GeneticMode::A, make_marker_variance());
    REQUIRE(gaussian.query<MarkerVarianceCap>() != nullptr);
    REQUIRE(gaussian.query<MixtureProportionCap>() == nullptr);
    REQUIRE(gaussian.query<MultiplierCap>() == nullptr);

    SpikeSlabGaussianPrior spike_slab(
        GeneticMode::A, make_marker_variance(), make_proportion_2());
    const auto* spike_slab_proportion
        = spike_slab.query<MixtureProportionCap>();
    REQUIRE(spike_slab.query<MarkerVarianceCap>() != nullptr);
    REQUIRE(spike_slab_proportion != nullptr);
    REQUIRE(spike_slab_proportion->proportion().size() == 1);
    REQUIRE(spike_slab.query<MultiplierCap>() == nullptr);

    ScaledMixtureGaussianPrior scaled_mixture(
        GeneticMode::A,
        make_marker_variance(),
        make_multiplier_2(),
        make_proportion_2());
    const auto* scaled_mixture_multiplier
        = scaled_mixture.query<MultiplierCap>();
    REQUIRE(scaled_mixture.query<MarkerVarianceCap>() != nullptr);
    REQUIRE(scaled_mixture.query<MixtureProportionCap>() != nullptr);
    REQUIRE(scaled_mixture_multiplier != nullptr);
    REQUIRE(scaled_mixture_multiplier->multiplier().size() == 1);

    JointMixtureGaussianPrior joint(
        std::array<GeneticMode, 2>{GeneticMode::A, GeneticMode::D},
        std::array<MarkerVariance, 2>{
            make_marker_variance(), make_marker_variance()},
        make_proportion_2());
    const auto* joint_variance = joint.query<MarkerVarianceCap>();
    const auto* joint_proportion = joint.query<MixtureProportionCap>();
    REQUIRE(joint_variance != nullptr);
    REQUIRE(joint_variance->variance().size() == 2);
    REQUIRE(joint_proportion != nullptr);
    REQUIRE(joint_proportion->proportion().size() == 1);
    REQUIRE(joint.query<MultiplierCap>() == nullptr);
}

TEST_CASE(
    "GeneticPrior::make_state builds capability-composed state",
    "[bayes_prior]")
{
    GaussianPrior gaussian(GeneticMode::A, make_marker_variance());
    constexpr Eigen::Index kSingleNumMarkers = 3;
    constexpr Eigen::Index kNumIndividuals = 5;

    auto gaussian_state
        = gaussian.make_state(kSingleNumMarkers, kNumIndividuals);
    auto& gaussian_variance = gaussian_state->require<VarianceStateCap>();
    REQUIRE(gaussian_variance.variance().size() == 1);
    REQUIRE(gaussian_variance.variance()[0].size() == 3);
    REQUIRE(gaussian_variance.variance()[0].isApprox(
        Eigen::VectorXd::Constant(3, 1.0)));
    REQUIRE(gaussian_state->query<ProportionStateCap>() == nullptr);

    SpikeSlabGaussianPrior spike_slab(
        GeneticMode::A, make_marker_variance(), make_proportion_2());
    auto spike_slab_state
        = spike_slab.make_state(kSingleNumMarkers, kNumIndividuals);
    auto& spike_slab_variance = spike_slab_state->require<VarianceStateCap>();
    auto& spike_slab_proportion
        = spike_slab_state->require<ProportionStateCap>();
    REQUIRE(spike_slab_variance.variance().size() == 1);
    REQUIRE(spike_slab_proportion.proportion().size() == 1);
    REQUIRE(spike_slab_proportion.proportion()[0].assignment.size() == 3);
    REQUIRE(spike_slab_proportion.proportion()[0].assignment.isZero());
    REQUIRE(spike_slab_proportion.proportion()[0].count.isApprox(
        Eigen::VectorXi{{3, 0}}));
    REQUIRE(spike_slab_proportion.proportion()[0].value.isApprox(
        Eigen::VectorXd{{0.9, 0.1}}));

    ScaledMixtureGaussianPrior scaled_mixture(
        GeneticMode::A,
        make_marker_variance(),
        make_multiplier_2(),
        make_proportion_2());
    auto scaled_mixture_state
        = scaled_mixture.make_state(kSingleNumMarkers, kNumIndividuals);
    auto& scaled_mixture_component
        = scaled_mixture_state->require<ComponentStateCap>();
    auto& scaled_mixture_proportion
        = scaled_mixture_state->require<ProportionStateCap>();
    REQUIRE(scaled_mixture_component.component().size() == 1);
    REQUIRE(scaled_mixture_component.component()[0].gebv.size() == 1);
    REQUIRE(
        scaled_mixture_component.component()[0].gebv[0].size()
        == kNumIndividuals);
    REQUIRE(scaled_mixture_component.component()[0].gebv_var.size() == 1);
    REQUIRE(scaled_mixture_proportion.proportion().size() == 1);
    REQUIRE(scaled_mixture_proportion.proportion()[0].assignment.size() == 3);
    auto& scaled_mixture_variance
        = scaled_mixture_state->require<VarianceStateCap>();
    REQUIRE(scaled_mixture_variance.variance().size() == 1);
    REQUIRE(scaled_mixture_variance.variance()[0].size() == kSingleNumMarkers);

    JointMixtureGaussianPrior joint(
        std::array<GeneticMode, 2>{GeneticMode::A, GeneticMode::D},
        std::array<MarkerVariance, 2>{
            make_marker_variance(), make_marker_variance()},
        make_proportion_2());

    auto joint_state = joint.make_state(4, kNumIndividuals);
    auto& joint_variance = joint_state->require<VarianceStateCap>();
    auto& joint_component = joint_state->require<ComponentStateCap>();
    auto& joint_proportion = joint_state->require<ProportionStateCap>();
    REQUIRE(joint_variance.variance().size() == 2);
    REQUIRE(joint_variance.variance()[0].size() == 4);
    REQUIRE(joint_variance.variance()[1].size() == 4);
    REQUIRE(joint_component.component().size() == 1);
    REQUIRE(joint_component.component()[0].gebv.size() == 1);
    REQUIRE(joint_component.component()[0].gebv[0].size() == kNumIndividuals);
    REQUIRE(joint_component.component()[0].gebv_var.size() == 1);
    REQUIRE(joint_proportion.proportion().size() == 1);
    REQUIRE(joint_proportion.proportion()[0].assignment.size() == 4);
}

TEST_CASE("BayesPrior::genetics range supports for-each", "[bayes_prior]")
{
    std::vector<std::unique_ptr<GeneticPrior>> genetics;
    genetics.push_back(make_gaussian(GeneticMode::A));
    genetics.push_back(make_spike_slab(GeneticMode::D));

    BayesPrior prior(
        make_random_prior(), std::move(genetics), make_residual_prior());

    std::vector<GeneticMode> visited;
    for (const auto& block : prior.genetics())
    {
        for (const auto m : block.modes())
        {
            visited.push_back(m);
        }
    }
    REQUIRE(
        visited == std::vector<GeneticMode>{GeneticMode::A, GeneticMode::D});
}

TEST_CASE(
    "BayesPrior::genetics empty for prior with no blocks",
    "[bayes_prior]")
{
    BayesPrior prior(make_random_prior(), {}, make_residual_prior());

    const auto range = prior.genetics();
    REQUIRE(range.empty());
    REQUIRE(range.size() == 0);
    REQUIRE(range.begin() == range.end());
}

TEST_CASE("BayesPrior rejects null genetic block", "[bayes_prior]")
{
    std::vector<std::unique_ptr<GeneticPrior>> genetics;
    genetics.push_back(nullptr);

    REQUIRE_THROWS_AS(
        BayesPrior(
            make_random_prior(), std::move(genetics), make_residual_prior()),
        GelexException);
}

TEST_CASE("BayesPrior rejects duplicate mode across blocks", "[bayes_prior]")
{
    std::vector<std::unique_ptr<GeneticPrior>> genetics;
    genetics.push_back(make_gaussian(GeneticMode::A));
    genetics.push_back(make_gaussian(GeneticMode::A));

    REQUIRE_THROWS_AS(
        BayesPrior(
            make_random_prior(), std::move(genetics), make_residual_prior()),
        GelexException);
}

TEST_CASE(
    "VarianceParameter rejects invalid initial_value",
    "[variance_parameter]")
{
    SECTION("initial_value <= 0")
    {
        REQUIRE_THROWS_AS(
            VarianceParameter(0.0, ScaledInvChiSqPrior{4.0, 1.0}),
            GelexException);
        REQUIRE_THROWS_AS(
            VarianceParameter(-1.0, ScaledInvChiSqPrior{4.0, 1.0}),
            GelexException);
    }
}

TEST_CASE(
    "VarianceParameter accepts positive initial_value",
    "[variance_parameter]")
{
    const VarianceParameter parameter{2.0, ScaledInvChiSqPrior{4.0, 1.0}};

    REQUIRE(parameter.initial_value() == 2.0);
}

TEST_CASE(
    "SimplexParameter rejects invalid initial_value",
    "[simplex_parameter]")
{
    SECTION("initial_value size < 2")
    {
        REQUIRE_THROWS_AS(
            SimplexParameter(
                Eigen::VectorXd{{1.0}},
                DirichletPrior{Eigen::VectorXd{{1.0, 1.0}}}),
            GelexException);
    }

    SECTION("initial_value entries must be positive")
    {
        REQUIRE_THROWS_AS(
            SimplexParameter(
                Eigen::VectorXd{{1.0, 0.0}},
                DirichletPrior{Eigen::VectorXd{{1.0, 1.0}}}),
            GelexException);
    }

    SECTION("initial_value entries must sum to one")
    {
        REQUIRE_THROWS_AS(
            SimplexParameter(
                Eigen::VectorXd{{0.7, 0.2}},
                DirichletPrior{Eigen::VectorXd{{1.0, 1.0}}}),
            GelexException);
    }

    SECTION("prior size must match initial_value")
    {
        REQUIRE_THROWS_AS(
            SimplexParameter(
                Eigen::VectorXd{{0.5, 0.5}},
                DirichletPrior{Eigen::VectorXd{{1.0, 1.0, 1.0}}}),
            GelexException);
    }
}

TEST_CASE("DirichletPrior rejects invalid concentration", "[dirichlet_prior]")
{
    SECTION("concentration entries must be positive")
    {
        REQUIRE_THROWS_AS(
            DirichletPrior(Eigen::VectorXd{{1.0, 0.0}}), GelexException);
    }

    SECTION("concentration size < 2")
    {
        REQUIRE_THROWS_AS(
            DirichletPrior(Eigen::VectorXd{{1.0}}), GelexException);
    }
}

TEST_CASE("ScaledInvChiSqPrior rejects invalid inputs", "[chisq_prior]")
{
    SECTION("degrees_of_freedom <= 0")
    {
        REQUIRE_THROWS_AS(ScaledInvChiSqPrior(0.0, 1.0), GelexException);
    }
    SECTION("scale <= 0 for proper prior")
    {
        REQUIRE_THROWS_AS(ScaledInvChiSqPrior(4.0, 0.0), GelexException);
    }
}

TEST_CASE(
    "ScaledInvChiSqPrior accepts flat prior sentinel {-2, 0}",
    "[chisq_prior]")
{
    SECTION("flat sentinel is accepted")
    {
        REQUIRE_NOTHROW(ScaledInvChiSqPrior(-2.0, 0.0));
    }
}

TEST_CASE(
    "ScaledInvChiSqPrior rejects near-flat but invalid prior combos",
    "[chisq_prior]")
{
    SECTION("{-2, 1}: df matches sentinel but scale != 0")
    {
        REQUIRE_THROWS_AS(ScaledInvChiSqPrior(-2.0, 1.0), GelexException);
    }
    SECTION("{4, 0}: scale matches sentinel but df != -2")
    {
        REQUIRE_THROWS_AS(ScaledInvChiSqPrior(4.0, 0.0), GelexException);
    }
    SECTION("{-1, 0}: negative df but not sentinel, scale == 0")
    {
        REQUIRE_THROWS_AS(ScaledInvChiSqPrior(-1.0, 0.0), GelexException);
    }
}
