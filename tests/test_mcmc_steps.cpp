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
#include <cmath>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/algo/infer/mcmc/chain.h"
#include "gelex/algo/infer/mcmc/step.h"
#include "gelex/algo/infer/mcmc/steps/genetic_step.h"
#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/gaussian_prior_state.h"
#include "gelex/bayes/genetic/parameters.h"
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

auto make_design(gelex::GeneticMode mode = gelex::GeneticMode::A)
    -> gelex::bayes::GeneticDesign
{
    Eigen::MatrixXd data{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}};
    auto mean = data.colwise().mean().transpose().eval();
    auto stddev = Eigen::VectorXd::Ones(data.cols());
    return gelex::bayes::GeneticDesign{
        mode,
        gelex::test::GenotypeBuilder::build(
            std::move(data), std::move(mean), std::move(stddev))};
}

auto make_proportion(Eigen::Index size = 2) -> gelex::bayes::MixtureProportion
{
    auto initial
        = Eigen::VectorXd::Constant(size, 1.0 / static_cast<double>(size));
    auto alpha = Eigen::VectorXd::Ones(size);
    return gelex::bayes::MixtureProportion{gelex::bayes::SimplexParameter{
        std::move(initial), gelex::bayes::DirichletPrior{std::move(alpha)}}};
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

TEST_CASE("Single shared Gaussian step updates fused state", "[mcmc]")
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
    auto& shared_state = std::get<gelex::bayes::SingleSharedGaussianState>(
        block.prior_state());

    std::mt19937_64 rng{123};
    gelex::mcmc::SingleSharedGaussianStep step{
        design, prior, block, residual, rng};
    step.step();

    REQUIRE(block.state().coeffs.allFinite());
    REQUIRE(std::isfinite(block.state().variance));
    REQUIRE(std::isfinite(shared_state.variance()));
    REQUIRE(shared_state.variance() > 0.0);
}

TEST_CASE("Single per-marker Gaussian step updates fused state", "[mcmc]")
{
    auto design = make_design();
    gelex::bayes::SingleGeneticPrior prior{
        gelex::bayes::SinglePerMarkerGaussianPrior{
            gelex::GeneticMode::A,
            gelex::bayes::PerMarkerVariance{make_variance(0.1)}}};
    gelex::bayes::SingleGeneticBlockState block{design, prior};
    block.state().coeffs = Eigen::VectorXd::Zero(2);
    gelex::bayes::ResidualState residual{
        .y_adj = Eigen::VectorXd{{1.0, -0.5, 0.25}}, .variance = 1.0};
    auto& per_marker_state
        = std::get<gelex::bayes::SinglePerMarkerGaussianState>(
            block.prior_state());

    std::mt19937_64 rng{123};
    gelex::mcmc::SinglePerMarkerGaussianStep step{
        design, prior, block, residual, rng};
    step.step();

    REQUIRE(block.state().coeffs.allFinite());
    REQUIRE(std::isfinite(block.state().variance));
    REQUIRE(per_marker_state.variance().allFinite());
    REQUIRE((per_marker_state.variance().array() > 0.0).all());
}

TEST_CASE("Single shared spike-slab step updates fused state", "[mcmc]")
{
    auto design = make_design();
    gelex::bayes::SingleGeneticPrior prior{
        gelex::bayes::SingleSharedSpikeSlabGaussianPrior{
            gelex::GeneticMode::A,
            gelex::bayes::SharedMarkerVariance{make_variance(0.1)},
            make_proportion()}};
    gelex::bayes::SingleGeneticBlockState block{design, prior};
    gelex::bayes::ResidualState residual{
        .y_adj = Eigen::VectorXd{{1.0, -0.5, 0.25}}, .variance = 1.0};
    auto& mixture_state
        = std::get<gelex::bayes::SingleSharedSpikeSlabGaussianState>(
            block.prior_state());

    std::mt19937_64 rng{123};
    gelex::mcmc::SingleSharedSpikeSlabStep step{
        design, prior, block, residual, rng};
    step.step();

    REQUIRE(block.state().coeffs.allFinite());
    REQUIRE(std::isfinite(mixture_state.variance()));
    REQUIRE(mixture_state.variance() > 0.0);
    REQUIRE(mixture_state.proportion().allFinite());
    REQUIRE(std::abs(mixture_state.proportion().sum() - 1.0) < 1e-12);
    REQUIRE((mixture_state.assignment().array() >= 0).all());
    REQUIRE((mixture_state.assignment().array() <= 1).all());
}

TEST_CASE("Single per-marker spike-slab step updates fused state", "[mcmc]")
{
    auto design = make_design();
    gelex::bayes::SingleGeneticPrior prior{
        gelex::bayes::SinglePerMarkerSpikeSlabGaussianPrior{
            gelex::GeneticMode::A,
            gelex::bayes::PerMarkerVariance{make_variance(0.1)},
            make_proportion()}};
    gelex::bayes::SingleGeneticBlockState block{design, prior};
    gelex::bayes::ResidualState residual{
        .y_adj = Eigen::VectorXd{{1.0, -0.5, 0.25}}, .variance = 1.0};
    auto& mixture_state
        = std::get<gelex::bayes::SinglePerMarkerSpikeSlabGaussianState>(
            block.prior_state());

    std::mt19937_64 rng{123};
    gelex::mcmc::SinglePerMarkerSpikeSlabStep step{
        design, prior, block, residual, rng};
    step.step();

    REQUIRE(block.state().coeffs.allFinite());
    REQUIRE(mixture_state.variance().allFinite());
    REQUIRE((mixture_state.variance().array() > 0.0).all());
    REQUIRE(mixture_state.proportion().allFinite());
    REQUIRE(std::abs(mixture_state.proportion().sum() - 1.0) < 1e-12);
    REQUIRE((mixture_state.assignment().array() >= 0).all());
    REQUIRE((mixture_state.assignment().array() <= 1).all());
}

TEST_CASE("Single scaled mixture step updates fused state", "[mcmc]")
{
    auto design = make_design();
    gelex::bayes::SingleGeneticPrior prior{
        gelex::bayes::SingleScaledMixtureGaussianPrior{
            gelex::GeneticMode::A,
            gelex::bayes::SharedMarkerVariance{make_variance(0.1)},
            Eigen::VectorXd{{0.0, 2.0}},
            make_proportion()}};
    gelex::bayes::SingleGeneticBlockState block{design, prior};
    gelex::bayes::ResidualState residual{
        .y_adj = Eigen::VectorXd{{1.0, -0.5, 0.25}}, .variance = 1.0};
    auto& scaled_state
        = std::get<gelex::bayes::SingleScaledMixtureGaussianState>(
            block.prior_state());

    std::mt19937_64 rng{123};
    gelex::mcmc::SingleScaledMixtureStep step{
        design, prior, block, residual, rng};
    step.step();

    REQUIRE(block.state().coeffs.allFinite());
    REQUIRE(std::isfinite(scaled_state.variance()));
    REQUIRE(scaled_state.variance() > 0.0);
    REQUIRE(scaled_state.component().gebv_var.allFinite());
    REQUIRE(scaled_state.proportion().allFinite());
    REQUIRE(std::abs(scaled_state.proportion().sum() - 1.0) < 1e-12);
    REQUIRE((scaled_state.assignment().array() >= 0).all());
    REQUIRE((scaled_state.assignment().array() <= 1).all());
}

TEST_CASE("Joint Gaussian mixture step updates fused state", "[mcmc]")
{
    auto additive = make_design(gelex::GeneticMode::A);
    auto dominance = make_design(gelex::GeneticMode::D);
    gelex::bayes::JointGeneticPrior prior{
        gelex::bayes::JointGaussianMixturePrior{
            gelex::bayes::JointSharedMarkerVariance{std::array{
                gelex::bayes::SharedMarkerVariance{make_variance(0.1)},
                gelex::bayes::SharedMarkerVariance{make_variance(0.2)}}},
            make_proportion(4)}};
    gelex::bayes::JointGeneticBlockState block{additive, dominance, prior};
    gelex::bayes::ResidualState residual{
        .y_adj = Eigen::VectorXd{{1.0, -0.5, 0.25}}, .variance = 1.0};
    auto& mixture_state = std::get<gelex::bayes::JointGaussianMixtureState>(
        block.prior_state());

    std::mt19937_64 rng{123};
    gelex::mcmc::JointGaussianMixtureStep step{
        additive, dominance, prior, block, residual, rng};
    step.step();

    REQUIRE(block.state(gelex::GeneticMode::A).coeffs.allFinite());
    REQUIRE(block.state(gelex::GeneticMode::D).coeffs.allFinite());
    REQUIRE(std::isfinite(mixture_state.variance(gelex::GeneticMode::A)));
    REQUIRE(std::isfinite(mixture_state.variance(gelex::GeneticMode::D)));
    REQUIRE(mixture_state.variance(gelex::GeneticMode::A) > 0.0);
    REQUIRE(mixture_state.variance(gelex::GeneticMode::D) > 0.0);
    REQUIRE(mixture_state.proportion().allFinite());
    REQUIRE(std::abs(mixture_state.proportion().sum() - 1.0) < 1e-12);
    REQUIRE((mixture_state.assignment().array() >= 0).all());
    REQUIRE((mixture_state.assignment().array() <= 3).all());
}
