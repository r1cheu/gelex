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
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/algo/infer/mcmc/chain.h"
#include "gelex/algo/infer/mcmc/step.h"
#include "gelex/algo/infer/mcmc/steps/genetic_variance.h"
#include "gelex/model/bayes/genetic_priors/gaussian.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/state.h"
#include "gelex/model/bayes/state_capabilities.h"
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

auto make_proportion() -> gelex::bayes::MixtureProportion
{
    return gelex::bayes::MixtureProportion{
        gelex::bayes::SimplexParameter{
            Eigen::VectorXd{{0.5, 0.5}},
            gelex::bayes::DirichletPrior{Eigen::VectorXd::Ones(2)}},
        gelex::bayes::UpdatePolicy::sampled};
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
    gelex::bayes::SingleSharedGaussianPrior prior{
        gelex::GeneticMode::A,
        gelex::bayes::SharedMarkerVariance{make_variance(0.1)}};
    gelex::bayes::SingleGeneticBlockState block{design, prior};
    block.state().coeffs = Eigen::VectorXd{{0.2, -0.3}};

    std::mt19937_64 rng{123};
    gelex::mcmc::SingleSharedGeneticVarianceStep step{prior, block, rng};
    step.step();

    const double variance
        = block.prior_state()
              .get<gelex::bayes::SingleSharedVarianceStateCap>()
              .variance();
    REQUIRE(std::isfinite(variance));
    REQUIRE(variance > 0.0);
}

TEST_CASE("Single scaled mixture variance step uses multiplier capability", "[mcmc]")
{
    auto design = make_design();
    gelex::bayes::SingleScaledMixtureGaussianPrior prior{
        gelex::GeneticMode::A,
        gelex::bayes::SharedMarkerVariance{make_variance(0.1)},
        Eigen::VectorXd{{0.0, 2.0}},
        make_proportion()};
    gelex::bayes::SingleGeneticBlockState block{design, prior};
    block.state().coeffs = Eigen::VectorXd{{0.2, -0.3}};
    auto& proportion
        = block.prior_state()
              .get<gelex::bayes::SingleProportionStateCap>()
              .proportion();
    proportion.assignment = Eigen::VectorXi{{1, 1}};
    proportion.count = Eigen::VectorXi{{0, 2}};

    std::mt19937_64 rng{123};
    gelex::mcmc::SingleSharedScaledMixtureGeneticVarianceStep step{
        prior, block, rng};
    step.step();

    const double variance
        = block.prior_state()
              .get<gelex::bayes::SingleSharedVarianceStateCap>()
              .variance();
    REQUIRE(std::isfinite(variance));
    REQUIRE(variance > 0.0);
}
