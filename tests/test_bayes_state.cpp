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
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/half_normal_prior.h"
#include "gelex/bayes/genetic/half_normal_prior_state.h"
#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"
#include "gelex/data/genotype.h"
#include "gelex/exception.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_effect_type.h"
#include "genotype_fixture.h"

namespace
{

auto make_variance(double value) -> gelex::bayes::VarianceParameter
{
    return gelex::bayes::VarianceParameter{
        value, gelex::bayes::ScaledInvChiSqPrior{4.0, 1.0}};
}

auto make_genotype(Eigen::MatrixXd data) -> gelex::Genotype
{
    auto mean = data.colwise().mean().transpose().eval();
    auto stddev = Eigen::VectorXd::Ones(data.cols());
    return gelex::test::GenotypeBuilder::build(
        std::move(data), std::move(mean), std::move(stddev));
}

auto make_model() -> gelex::BayesModel
{
    std::vector<gelex::bayes::GeneticDesign> genetics;
    genetics.emplace_back(
        gelex::GeneticMode::A,
        make_genotype(Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}}));
    genetics.emplace_back(
        gelex::GeneticMode::D,
        make_genotype(Eigen::MatrixXd{{1.0, 0.0}, {0.0, 1.0}, {1.0, 2.0}}));

    return gelex::BayesModel{
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        gelex::FixedDesign::make(3),
        {},
        std::move(genetics)};
}

auto make_prior() -> gelex::bayes::BayesPrior
{
    std::vector<gelex::bayes::GeneticPrior> genetics;
    genetics.emplace_back(
        gelex::bayes::SingleGeneticPrior{
            gelex::bayes::SingleSharedGaussianPrior{
                gelex::GeneticMode::A,
                gelex::bayes::SharedMarkerVariance{make_variance(0.1)}}});
    genetics.emplace_back(
        gelex::bayes::SingleGeneticPrior{
            gelex::bayes::SinglePerMarkerGaussianPrior{
                gelex::GeneticMode::D,
                gelex::bayes::PerMarkerVariance{make_variance(0.2)}}});

    return gelex::bayes::BayesPrior{
        gelex::bayes::RandomPrior{make_variance(0.3)},
        std::move(genetics),
        gelex::bayes::ResidualPrior{make_variance(0.4)}};
}

}  // namespace

TEST_CASE("BayesState creates single genetic blocks", "[bayes_state]")
{
    auto model = make_model();
    auto prior = make_prior();
    gelex::BayesState state(model, prior);

    REQUIRE(state.genetics().size() == 2);
    REQUIRE(
        std::holds_alternative<gelex::bayes::SingleGeneticBlockState>(
            state.genetics()[0]));
    REQUIRE(
        std::holds_alternative<gelex::bayes::SingleGeneticBlockState>(
            state.genetics()[1]));
    REQUIRE(state.residual().variance == 0.4);
}

TEST_CASE("BayesState creates joint half normal mixture block", "[bayes_state]")
{
    auto model = make_model();
    std::vector<gelex::bayes::GeneticPrior> genetics;
    genetics.emplace_back(
        gelex::bayes::JointGeneticPrior{
            gelex::bayes::JointHalfNormalMixturePrior{
                gelex::bayes::JointSharedMarkerVariance{std::array{
                    gelex::bayes::SharedMarkerVariance{make_variance(0.1)},
                    gelex::bayes::SharedMarkerVariance{make_variance(0.2)}}},
                gelex::bayes::MixtureProportion{
                    Eigen::VectorXd{{0.7, 0.1, 0.1, 0.1}}},
                gelex::bayes::ProbabilityParameter{
                    0.6, gelex::bayes::BetaPrior{1.0, 1.0}}}});
    gelex::bayes::BayesPrior prior{
        gelex::bayes::RandomPrior{make_variance(0.3)},
        std::move(genetics),
        gelex::bayes::ResidualPrior{make_variance(0.4)}};

    gelex::BayesState state(model, prior);
    auto& block
        = std::get<gelex::bayes::JointGeneticBlockState>(state.genetics()[0]);
    auto& prior_state = std::get<gelex::bayes::JointHalfNormalMixtureState>(
        block.prior_state());

    REQUIRE(state.genetics().size() == 1);
    REQUIRE(prior_state.variance(gelex::GeneticMode::A) == 0.1);
    REQUIRE(prior_state.variance(gelex::GeneticMode::D) == 0.2);
    REQUIRE(prior_state.dominance_sign().positive_probability == 0.6);
    REQUIRE(prior_state.dominance_sign().sign.size() == 2);
}

TEST_CASE("BayesState rejects missing genetic designs", "[bayes_state]")
{
    auto model = make_model();
    std::vector<gelex::bayes::GeneticPrior> genetics;
    genetics.emplace_back(
        gelex::bayes::SingleGeneticPrior{
            gelex::bayes::SingleSharedGaussianPrior{
                gelex::GeneticMode::A,
                gelex::bayes::SharedMarkerVariance{make_variance(0.1)}}});
    genetics.emplace_back(
        gelex::bayes::JointGeneticPrior{gelex::bayes::JointGaussianMixturePrior{
            gelex::bayes::JointSharedMarkerVariance{std::array{
                gelex::bayes::SharedMarkerVariance{make_variance(0.1)},
                gelex::bayes::SharedMarkerVariance{make_variance(0.2)}}},
            gelex::bayes::MixtureProportion{
                Eigen::VectorXd{{0.25, 0.25, 0.25, 0.25}}}}});

    REQUIRE_THROWS_AS(
        gelex::bayes::BayesPrior(
            gelex::bayes::RandomPrior{make_variance(0.3)},
            std::move(genetics),
            gelex::bayes::ResidualPrior{make_variance(0.4)}),
        gelex::GelexException);
    static_cast<void>(model);
}
