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

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/exception.h"
#include "gelex/types/constrained_value.h"
#include "gelex/types/constrained_vector.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/genetic_priors/gaussian.h"
#include "gelex/model/bayes/genetic_priors/joint_mixture_gaussian.h"
#include "gelex/model/bayes/genetic_priors/mixture_gaussian.h"
#include "gelex/model/bayes/legacy_prior.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/prior_specs.h"
#include "gelex/types/genetic_effect_type.h"

using gelex::GelexException;
using gelex::GeneticMode;
using gelex::bayes::BayesPrior;
using gelex::bayes::GaussianPrior;
using gelex::bayes::GeneticPrior;
using gelex::bayes::JointMixtureGaussianPrior;
using gelex::bayes::MarkerVarianceScope;
using gelex::bayes::MarkerVarianceSpec;
using gelex::PositiveScalar;
using gelex::PositiveVector;
using gelex::bayes::ProportionSpec;
using gelex::bayes::ProportionUpdate;
using gelex::bayes::ScaledInvChiSqPrior;
using gelex::Simplex;
using gelex::bayes::SpikeSlabGaussianPrior;
using gelex::bayes::VarianceSpec;

namespace
{

auto make_variance(double init = 1.0) -> VarianceSpec
{
    return VarianceSpec(init, ScaledInvChiSqPrior{4.0, 1.0});
}

auto make_marker_variance() -> MarkerVarianceSpec
{
    return MarkerVarianceSpec{
        MarkerVarianceScope::per_marker,
        make_variance(),
    };
}

auto make_proportion_2() -> ProportionSpec
{
    return ProportionSpec{
        Simplex<double>{{0.9, 0.1}},
        PositiveVector<double>{{1.0, 1.0}},
        ProportionUpdate::fixed,
    };
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
        std::array<MarkerVarianceSpec, 2>{
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
        make_variance(2.0), std::move(genetics), make_variance(5.0));

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

TEST_CASE("BayesPrior::genetics range supports for-each", "[bayes_prior]")
{
    std::vector<std::unique_ptr<GeneticPrior>> genetics;
    genetics.push_back(make_gaussian(GeneticMode::A));
    genetics.push_back(make_spike_slab(GeneticMode::D));

    BayesPrior prior(make_variance(), std::move(genetics), make_variance());

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
    BayesPrior prior(make_variance(), {}, make_variance());

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
        BayesPrior(make_variance(), std::move(genetics), make_variance()),
        GelexException);
}

TEST_CASE("BayesPrior rejects duplicate mode across blocks", "[bayes_prior]")
{
    std::vector<std::unique_ptr<GeneticPrior>> genetics;
    genetics.push_back(make_gaussian(GeneticMode::A));
    genetics.push_back(make_gaussian(GeneticMode::A));

    REQUIRE_THROWS_AS(
        BayesPrior(make_variance(), std::move(genetics), make_variance()),
        GelexException);
}

TEST_CASE("VarianceSpec rejects invalid initial_value", "[variance_spec]")
{
    SECTION("initial_value <= 0")
    {
        REQUIRE_THROWS_AS(
            VarianceSpec(0.0, ScaledInvChiSqPrior{4.0, 1.0}), GelexException);
        REQUIRE_THROWS_AS(
            VarianceSpec(-1.0, ScaledInvChiSqPrior{4.0, 1.0}), GelexException);
    }
}

TEST_CASE("VarianceSpec accepts constrained initial_value", "[variance_spec]")
{
    const VarianceSpec spec{
        PositiveScalar<double>{2.0}, ScaledInvChiSqPrior{4.0, 1.0}};

    REQUIRE(spec.initial_value() == 2.0);
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
