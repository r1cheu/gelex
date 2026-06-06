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
#include <utility>
#include <variant>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "gelex/exception.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/genetic/gaussian_prior_state.h"
#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/half_normal_prior.h"
#include "gelex/bayes/genetic/half_normal_prior_state.h"
#include "gelex/bayes/prior.h"
#include "gelex/types/genetic_effect_type.h"

namespace
{

auto make_variance(double value) -> gelex::bayes::VarianceParameter
{
    return gelex::bayes::VarianceParameter{
        value, gelex::bayes::ScaledInvChiSqPrior{4.0, 1.0}};
}

auto make_random_prior() -> gelex::bayes::RandomPrior
{
    return gelex::bayes::RandomPrior{make_variance(0.2)};
}

auto make_residual_prior() -> gelex::bayes::ResidualPrior
{
    return gelex::bayes::ResidualPrior{make_variance(0.8)};
}

auto make_shared_prior(gelex::GeneticMode mode) -> gelex::bayes::GeneticPrior
{
    return gelex::bayes::SingleGeneticPrior{
        gelex::bayes::SingleSharedGaussianPrior{
            mode, gelex::bayes::SharedMarkerVariance{make_variance(0.1)}}};
}

}  // namespace

TEST_CASE("BayesPrior owns single and joint genetic blocks", "[bayes_prior]")
{
    std::vector<gelex::bayes::GeneticPrior> genetics;
    genetics.emplace_back(make_shared_prior(gelex::GeneticMode::A));

    gelex::bayes::BayesPrior prior{
        make_random_prior(), std::move(genetics), make_residual_prior()};

    REQUIRE(prior.genetics().size() == 1);
    const auto& genetic
        = std::get<gelex::bayes::SingleGeneticPrior>(prior.genetics()[0]);
    REQUIRE(gelex::bayes::mode(genetic) == gelex::GeneticMode::A);
}

TEST_CASE("Single shared genetic prior creates concrete state", "[bayes_prior]")
{
    gelex::bayes::SingleSharedGaussianPrior prior{
        gelex::GeneticMode::A,
        gelex::bayes::SharedMarkerVariance{make_variance(0.1)}};
    auto state = prior.make_state(3, 2);

    REQUIRE(state.variance() == 0.1);
}

TEST_CASE("Probability parameter validates beta prior", "[bayes_prior]")
{
    REQUIRE_NOTHROW(
        (gelex::bayes::ProbabilityParameter{
            0.5, gelex::bayes::BetaPrior{1.0, 1.0}}));
    REQUIRE_THROWS_AS(
        (gelex::bayes::BetaPrior{0.0, 1.0}), gelex::GelexException);
    REQUIRE_THROWS_AS(
        (gelex::bayes::BetaPrior{1.0, 0.0}), gelex::GelexException);
    REQUIRE_THROWS_AS(
        (gelex::bayes::ProbabilityParameter{
            0.0, gelex::bayes::BetaPrior{1.0, 1.0}}),
        gelex::GelexException);
    REQUIRE_THROWS_AS(
        (gelex::bayes::ProbabilityParameter{
            1.0, gelex::bayes::BetaPrior{1.0, 1.0}}),
        gelex::GelexException);
}

TEST_CASE("Joint half normal mixture prior creates concrete state", "[bayes_prior]")
{
    gelex::bayes::JointHalfNormalMixturePrior prior{
        gelex::bayes::JointSharedMarkerVariance{std::array{
            gelex::bayes::SharedMarkerVariance{make_variance(0.1)},
            gelex::bayes::SharedMarkerVariance{make_variance(0.2)}}},
        gelex::bayes::MixtureProportion{
            Eigen::VectorXd{{0.7, 0.1, 0.1, 0.1}}},
        gelex::bayes::ProbabilityParameter{
            0.6, gelex::bayes::BetaPrior{1.0, 1.0}}};
    auto state = prior.make_state(3, 2);

    REQUIRE(state.variance(gelex::GeneticMode::A) == 0.1);
    REQUIRE(state.variance(gelex::GeneticMode::D) == 0.2);
    REQUIRE(state.proportion().isApprox(
        Eigen::VectorXd{{0.7, 0.1, 0.1, 0.1}}));
    REQUIRE((state.assignment().values().array() == 0).all());
    REQUIRE(state.dominance_sign().positive_probability == 0.6);
    REQUIRE((state.dominance_sign().sign.values().array() == 0).all());
}

TEST_CASE(
    "Joint half normal mixture prior rejects non-four-class proportion",
    "[bayes_prior]")
{
    REQUIRE_THROWS_AS(
        (gelex::bayes::JointHalfNormalMixturePrior{
            gelex::bayes::JointSharedMarkerVariance{std::array{
                gelex::bayes::SharedMarkerVariance{make_variance(0.1)},
                gelex::bayes::SharedMarkerVariance{make_variance(0.2)}}},
            gelex::bayes::MixtureProportion{Eigen::VectorXd{{0.5, 0.5}}},
            gelex::bayes::ProbabilityParameter{
                0.6, gelex::bayes::BetaPrior{1.0, 1.0}}}),
        gelex::GelexException);
}

TEST_CASE("BayesPrior rejects duplicate genetic modes", "[bayes_prior]")
{
    std::vector<gelex::bayes::GeneticPrior> genetics;
    genetics.emplace_back(make_shared_prior(gelex::GeneticMode::A));
    genetics.emplace_back(make_shared_prior(gelex::GeneticMode::A));

    REQUIRE_THROWS_AS(
        gelex::bayes::BayesPrior(
            make_random_prior(), std::move(genetics), make_residual_prior()),
        gelex::GelexException);
}
