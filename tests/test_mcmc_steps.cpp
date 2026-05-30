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

#include <cmath>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/algo/infer/mcmc/chain.h"
#include "gelex/algo/infer/mcmc/step.h"
#include "gelex/algo/infer/mcmc/steps/genetic_coefficient.h"
#include "gelex/algo/infer/mcmc/steps/genetic_proportion.h"
#include "gelex/algo/infer/mcmc/steps/genetic_variance.h"
#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/genetic/gaussian_prior_state.h"
#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/state.h"
#include "gelex/types/genetic_effect_type.h"
#include "genotype_fixture.h"

namespace
{

auto make_variance(double value) -> gelex::bayes::VarianceParameter
{
    return gelex::bayes::VarianceParameter{
        value, gelex::bayes::ScaledInvChiSqPrior{4.0, 1.0}};
}

auto make_design() -> gelex::bayes::GeneticDesign
{
    Eigen::MatrixXd data{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}};
    auto mean = data.colwise().mean().transpose().eval();
    auto stddev = Eigen::VectorXd::Ones(data.cols());
    return gelex::bayes::GeneticDesign{
        gelex::GeneticMode::A,
        gelex::test::GenotypeBuilder::build(
            std::move(data), std::move(mean), std::move(stddev))};
}

auto make_proportion() -> gelex::bayes::SampledProportion
{
    return gelex::bayes::SampledProportion{gelex::bayes::SimplexParameter{
        Eigen::VectorXd{{0.5, 0.5}},
        gelex::bayes::DirichletPrior{Eigen::VectorXd::Ones(2)}}};
}

class CountingStep final : public gelex::mcmc::Step
{
   public:
    explicit CountingStep(int& count) : count_(count) {}

    auto step() -> void override { ++count_; }

   private:
    int& count_;
};

}  // namespace

TEST_CASE("Runtime Chain runs heterogeneous steps in order", "[mcmc][chain]")
{
    int first = 0;
    int second = 0;
    std::vector<std::unique_ptr<gelex::mcmc::Step>> steps;
    steps.push_back(std::make_unique<CountingStep>(first));
    steps.push_back(std::make_unique<CountingStep>(second));

    gelex::mcmc::Chain chain{std::move(steps)};
    chain.step();

    REQUIRE(first == 1);
    REQUIRE(second == 1);
}

TEST_CASE("Single shared variance step samples finite variance", "[mcmc]")
{
    auto design = make_design();
    gelex::bayes::SingleGeneticPrior prior{
        gelex::bayes::SingleSharedGaussianPrior{
            gelex::GeneticMode::A,
            gelex::bayes::SharedMarkerVariance{make_variance(0.1)}}};
    gelex::bayes::SingleGeneticBlockState block{design, prior};
    block.state().coeffs = Eigen::VectorXd{{0.2, -0.3}};
    auto& shared_state = std::get<gelex::bayes::SingleSharedGaussianState>(
        block.prior_state());

    std::mt19937_64 rng{123};
    gelex::mcmc::SingleSharedVarStep step{block, prior, rng};
    step.step();

    const double variance = shared_state.variance();
    REQUIRE(std::isfinite(variance));
    REQUIRE(variance > 0.0);
}

TEST_CASE(
    "Single scaled mixture variance step uses multiplier capability",
    "[mcmc]")
{
    auto design = make_design();
    gelex::bayes::SingleGeneticPrior prior{
        gelex::bayes::SingleSampledScaledMixtureGaussianPrior{
            gelex::GeneticMode::A,
            gelex::bayes::SharedMarkerVariance{make_variance(0.1)},
            Eigen::VectorXd{{0.0, 2.0}},
            make_proportion()}};
    gelex::bayes::SingleGeneticBlockState block{design, prior};
    block.state().coeffs = Eigen::VectorXd{{0.2, -0.3}};
    auto& scaled_state
        = std::get<gelex::bayes::SingleSampledScaledMixtureGaussianState>(
            block.prior_state());
    scaled_state.assignment().assignment = Eigen::VectorXi{{1, 1}};
    scaled_state.assignment().count = Eigen::VectorXi{{0, 2}};

    std::mt19937_64 rng{123};
    gelex::mcmc::SingleSampledScaledMixtureVarStep step{block, prior, rng};
    step.step();

    const double variance = scaled_state.variance();
    REQUIRE(std::isfinite(variance));
    REQUIRE(variance > 0.0);
}

TEST_CASE("Single proportion step accepts prior state variant", "[mcmc]")
{
    auto design = make_design();
    gelex::bayes::SingleGeneticPrior prior{
        gelex::bayes::SingleSampledSharedSpikeSlabGaussianPrior{
            gelex::GeneticMode::A,
            gelex::bayes::SharedMarkerVariance{make_variance(0.1)},
            make_proportion()}};
    gelex::bayes::SingleGeneticBlockState block{design, prior};
    auto& proportion_state
        = std::get<gelex::bayes::SingleSampledSharedSpikeSlabGaussianState>(
              block.prior_state())
              .proportion();
    proportion_state.assignment.count = Eigen::VectorXi{{1, 1}};

    std::mt19937_64 rng{123};
    gelex::mcmc::SingleProportionStep step{block, prior, rng};
    step.step();

    REQUIRE(proportion_state.value.allFinite());
    REQUIRE(std::abs(proportion_state.value.sum() - 1.0) < 1e-12);
}

TEST_CASE("Single coefficient step accepts prior state variant", "[mcmc]")
{
    auto design = make_design();
    gelex::bayes::SingleGeneticPrior prior{
        gelex::bayes::SingleSharedGaussianPrior{
            gelex::GeneticMode::A,
            gelex::bayes::SharedMarkerVariance{make_variance(0.1)}}};
    gelex::bayes::SingleGeneticBlockState block{design, prior};
    block.state().coeffs = Eigen::VectorXd::Zero(2);
    gelex::bayes::ResidualState residual{
        .y_adj = Eigen::VectorXd{{1.0, -0.5, 0.25}}, .variance = 1.0};

    std::mt19937_64 rng{123};
    gelex::mcmc::SingleSharedGaussianCoeffStep step{
        design, block, prior, residual, rng};
    step.step();

    REQUIRE(block.state().coeffs.allFinite());
    REQUIRE(std::isfinite(block.state().variance));
}
