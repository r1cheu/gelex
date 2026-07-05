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
#include <cstddef>
#include <filesystem>
#include <random>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "file_fixture.h"
#include "gelex/algo/mcmc/chain.h"
#include "gelex/algo/mcmc/solver.h"
#include "gelex/algo/mcmc/steps/joint_genetic_step.h"
#include "gelex/algo/mcmc/steps/random.h"
#include "gelex/algo/mcmc/steps/single_genetic_step.h"
#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/gaussian_prior_state.h"
#include "gelex/bayes/genetic/half_normal_prior.h"
#include "gelex/bayes/genetic/half_normal_prior_state.h"
#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/infra/stats/result.h"
#include "gelex/io/binary_reader.h"
#include "gelex/io/mcmc.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_mode.h"
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

auto make_active_only_design(gelex::GeneticMode mode = gelex::GeneticMode::A)
    -> gelex::bayes::GeneticDesign
{
    Eigen::MatrixXd data{{0.0}, {1.0}, {2.0}};
    auto mean = data.colwise().mean().transpose().eval();
    auto stddev = Eigen::VectorXd::Ones(data.cols());
    return gelex::bayes::GeneticDesign{
        mode,
        gelex::test::GenotypeBuilder::build(
            std::move(data), std::move(mean), std::move(stddev))};
}

auto make_design_with_monomorphic_second_marker(
    gelex::GeneticMode mode = gelex::GeneticMode::A)
    -> gelex::bayes::GeneticDesign
{
    Eigen::MatrixXd data{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}};
    auto mean = data.colwise().mean().transpose().eval();
    auto stddev = Eigen::VectorXd::Ones(data.cols());
    return gelex::bayes::GeneticDesign{
        mode,
        gelex::test::GenotypeBuilder::build(
            std::move(data), std::move(mean), std::move(stddev), {1})};
}

auto make_proportion(Eigen::Index size = 2) -> gelex::bayes::MixtureProportion
{
    auto initial
        = Eigen::VectorXd::Constant(size, 1.0 / static_cast<double>(size));
    auto alpha = Eigen::VectorXd::Ones(size);
    return gelex::bayes::MixtureProportion{gelex::bayes::SimplexParameter{
        std::move(initial), gelex::bayes::DirichletPrior{std::move(alpha)}}};
}

}  // namespace

TEST_CASE("Random step updates coefficients and variance", "[mcmc]")
{
    std::vector<gelex::bayes::RandomDesign> designs;
    designs.emplace_back(
        "pen",
        std::vector<std::string>{"a", "b"},
        Eigen::MatrixXd{{1.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}});
    gelex::bayes::RandomPrior prior{make_variance(0.5)};
    std::vector<gelex::bayes::RandomState> states;
    states.emplace_back(designs.front(), prior);
    gelex::bayes::ResidualState residual{
        .y_adj = Eigen::VectorXd{{1.0, -0.5, 0.25}}, .variance = 1.0};

    std::mt19937_64 rng{123};
    gelex::mcmc::RandomStep step{prior, designs, states, residual, rng};
    step.step();

    REQUIRE(states.front().coeffs.allFinite());
    REQUIRE(std::isfinite(states.front().variance));
    REQUIRE(states.front().variance > 0.0);
    REQUIRE(residual.y_adj.allFinite());
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
    REQUIRE((mixture_state.assignment().values().array() >= 0).all());
    REQUIRE((mixture_state.assignment().values().array() <= 1).all());
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
    REQUIRE((mixture_state.assignment().values().array() >= 0).all());
    REQUIRE((mixture_state.assignment().values().array() <= 1).all());
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
    REQUIRE((scaled_state.assignment().values().array() >= 0).all());
    REQUIRE((scaled_state.assignment().values().array() <= 1).all());
}

TEST_CASE(
    "Single scaled mixture step skips monomorphic markers completely",
    "[mcmc]")
{
    auto with_mono = make_design_with_monomorphic_second_marker();
    auto active_only = make_active_only_design();
    gelex::bayes::SingleGeneticPrior with_mono_prior{
        gelex::bayes::SingleScaledMixtureGaussianPrior{
            gelex::GeneticMode::A,
            gelex::bayes::SharedMarkerVariance{make_variance(0.1)},
            Eigen::VectorXd{{0.0, 2.0}},
            make_proportion()}};
    gelex::bayes::SingleGeneticPrior active_only_prior{
        gelex::bayes::SingleScaledMixtureGaussianPrior{
            gelex::GeneticMode::A,
            gelex::bayes::SharedMarkerVariance{make_variance(0.1)},
            Eigen::VectorXd{{0.0, 2.0}},
            make_proportion()}};
    gelex::bayes::SingleGeneticBlockState with_mono_block{
        with_mono, with_mono_prior};
    gelex::bayes::SingleGeneticBlockState active_only_block{
        active_only, active_only_prior};
    with_mono_block.state().coeffs = Eigen::VectorXd{{0.0, 50.0}};
    active_only_block.state().coeffs = Eigen::VectorXd{{0.0}};
    gelex::bayes::ResidualState with_mono_residual{
        .y_adj = Eigen::VectorXd{{1.0, -0.5, 0.25}}, .variance = 1.0};
    gelex::bayes::ResidualState active_only_residual{
        .y_adj = Eigen::VectorXd{{1.0, -0.5, 0.25}}, .variance = 1.0};
    auto& with_mono_state
        = std::get<gelex::bayes::SingleScaledMixtureGaussianState>(
            with_mono_block.prior_state());
    auto& active_only_state
        = std::get<gelex::bayes::SingleScaledMixtureGaussianState>(
            active_only_block.prior_state());
    with_mono_state.assignment() = Eigen::VectorXi{{0, 1}};
    active_only_state.assignment() = Eigen::VectorXi{{0}};

    std::mt19937_64 with_mono_rng{123};
    std::mt19937_64 active_only_rng{123};
    gelex::mcmc::SingleScaledMixtureStep with_mono_step{
        with_mono,
        with_mono_prior,
        with_mono_block,
        with_mono_residual,
        with_mono_rng};
    gelex::mcmc::SingleScaledMixtureStep active_only_step{
        active_only,
        active_only_prior,
        active_only_block,
        active_only_residual,
        active_only_rng};
    with_mono_step.step();
    active_only_step.step();

    REQUIRE(
        with_mono_block.state().coeffs(0)
        == active_only_block.state().coeffs(0));
    REQUIRE(with_mono_block.state().coeffs(1) == 50.0);
    REQUIRE(
        with_mono_state.assignment()(0) == active_only_state.assignment()(0));
    REQUIRE(with_mono_state.assignment()(1) == 1);
    REQUIRE(with_mono_state.variance() == active_only_state.variance());
    REQUIRE(
        with_mono_state.proportion().isApprox(active_only_state.proportion()));
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
    REQUIRE((mixture_state.assignment().values().array() >= 0).all());
    REQUIRE((mixture_state.assignment().values().array() <= 3).all());
}

TEST_CASE("Joint half normal mixture step updates fused state", "[mcmc]")
{
    auto additive = make_design(gelex::GeneticMode::A);
    auto dominance = make_design(gelex::GeneticMode::D);
    auto initial = Eigen::VectorXd{{0.001, 0.001, 0.997, 0.001}};
    auto alpha = Eigen::VectorXd::Ones(4);
    gelex::bayes::JointGeneticPrior prior{
        gelex::bayes::JointHalfNormalMixturePrior{
            gelex::bayes::JointSharedMarkerVariance{std::array{
                gelex::bayes::SharedMarkerVariance{make_variance(0.1)},
                gelex::bayes::SharedMarkerVariance{make_variance(0.2)}}},
            gelex::bayes::MixtureProportion{gelex::bayes::SimplexParameter{
                std::move(initial),
                gelex::bayes::DirichletPrior{std::move(alpha)}}},
            gelex::bayes::ProbabilityParameter{
                0.8, gelex::bayes::BetaPrior{1.0, 1.0}}}};
    gelex::bayes::JointGeneticBlockState block{additive, dominance, prior};
    gelex::bayes::ResidualState residual{
        .y_adj = Eigen::VectorXd{{1.0, -0.5, 0.25}}, .variance = 1.0};
    auto& mixture_state = std::get<gelex::bayes::JointHalfNormalMixtureState>(
        block.prior_state());

    std::mt19937_64 rng{123};
    gelex::mcmc::JointHalfNormalMixtureStep step{
        additive, dominance, prior, block, residual, rng};
    step.step();

    const auto& additive_coeffs = block.state(gelex::GeneticMode::A).coeffs;
    const auto& dominance_coeffs = block.state(gelex::GeneticMode::D).coeffs;
    REQUIRE(additive_coeffs.allFinite());
    REQUIRE(dominance_coeffs.allFinite());
    REQUIRE(std::isfinite(mixture_state.variance(gelex::GeneticMode::A)));
    REQUIRE(std::isfinite(mixture_state.variance(gelex::GeneticMode::D)));
    REQUIRE(mixture_state.variance(gelex::GeneticMode::A) > 0.0);
    REQUIRE(mixture_state.variance(gelex::GeneticMode::D) > 0.0);
    REQUIRE(mixture_state.proportion().allFinite());
    REQUIRE(std::abs(mixture_state.proportion().sum() - 1.0) < 1e-12);
    REQUIRE((mixture_state.assignment().values().array() >= 0).all());
    REQUIRE((mixture_state.assignment().values().array() <= 3).all());
    REQUIRE((mixture_state.dominance_sign().sign.values().array() >= 0).all());
    REQUIRE((mixture_state.dominance_sign().sign.values().array() <= 1).all());
    REQUIRE(std::isfinite(mixture_state.dominance_sign().positive_probability));
    REQUIRE(mixture_state.dominance_sign().positive_probability > 0.0);
    REQUIRE(mixture_state.dominance_sign().positive_probability < 1.0);

    bool has_dominance_active = false;
    for (Eigen::Index i = 0; i < mixture_state.assignment().size(); ++i)
    {
        const int component = mixture_state.assignment()(i);
        if (component == 0 || component == 2)
        {
            REQUIRE(additive_coeffs(i) == 0.0);
        }
        if (component == 0 || component == 1)
        {
            REQUIRE(dominance_coeffs(i) == 0.0);
        }
        if (component == 2 || component == 3)
        {
            has_dominance_active = true;
            if (mixture_state.dominance_sign().sign(i) == 1)
            {
                REQUIRE(dominance_coeffs(i) > 0.0);
            }
            else
            {
                REQUIRE(dominance_coeffs(i) < 0.0);
            }
        }
    }
    REQUIRE(has_dominance_active);
}

TEST_CASE(
    "Chain::make maps single genetic priors to MCMC steps",
    "[mcmc][chain]")
{
    SECTION("shared Gaussian")
    {
        std::vector<gelex::bayes::GeneticDesign> genetics;
        genetics.push_back(make_design());
        gelex::BayesModel model{
            Eigen::VectorXd{{1.0, -0.5, 0.25}},
            gelex::FixedDesign::make(3),
            {},
            std::move(genetics)};
        std::vector<gelex::bayes::GeneticPrior> priors;
        priors.emplace_back(
            gelex::bayes::SingleGeneticPrior{
                gelex::bayes::SingleSharedGaussianPrior{
                    gelex::GeneticMode::A,
                    gelex::bayes::SharedMarkerVariance{make_variance(0.1)}}});
        gelex::bayes::BayesPrior prior{
            gelex::bayes::RandomPrior{make_variance(0.3)},
            std::move(priors),
            gelex::bayes::ResidualPrior{make_variance(0.4)}};
        gelex::BayesState state(model, prior);

        std::mt19937_64 rng{123};
        auto chain = gelex::mcmc::Chain::make(model, prior, state, rng);
        chain.step();

        const auto& block = std::get<gelex::bayes::SingleGeneticBlockState>(
            state.genetics()[0]);
        REQUIRE(block.state().coeffs.allFinite());
        REQUIRE(std::isfinite(block.state().heritability));
    }

    SECTION("per-marker Gaussian")
    {
        std::vector<gelex::bayes::GeneticDesign> genetics;
        genetics.push_back(make_design());
        gelex::BayesModel model{
            Eigen::VectorXd{{1.0, -0.5, 0.25}},
            gelex::FixedDesign::make(3),
            {},
            std::move(genetics)};
        std::vector<gelex::bayes::GeneticPrior> priors;
        priors.emplace_back(
            gelex::bayes::SingleGeneticPrior{
                gelex::bayes::SinglePerMarkerGaussianPrior{
                    gelex::GeneticMode::A,
                    gelex::bayes::PerMarkerVariance{make_variance(0.1)}}});
        gelex::bayes::BayesPrior prior{
            gelex::bayes::RandomPrior{make_variance(0.3)},
            std::move(priors),
            gelex::bayes::ResidualPrior{make_variance(0.4)}};
        gelex::BayesState state(model, prior);

        std::mt19937_64 rng{123};
        auto chain = gelex::mcmc::Chain::make(model, prior, state, rng);
        chain.step();

        const auto& block = std::get<gelex::bayes::SingleGeneticBlockState>(
            state.genetics()[0]);
        REQUIRE(block.state().coeffs.allFinite());
        REQUIRE(std::isfinite(block.state().heritability));
    }

    SECTION("shared spike-slab")
    {
        std::vector<gelex::bayes::GeneticDesign> genetics;
        genetics.push_back(make_design());
        gelex::BayesModel model{
            Eigen::VectorXd{{1.0, -0.5, 0.25}},
            gelex::FixedDesign::make(3),
            {},
            std::move(genetics)};
        std::vector<gelex::bayes::GeneticPrior> priors;
        priors.emplace_back(
            gelex::bayes::SingleGeneticPrior{
                gelex::bayes::SingleSharedSpikeSlabGaussianPrior{
                    gelex::GeneticMode::A,
                    gelex::bayes::SharedMarkerVariance{make_variance(0.1)},
                    make_proportion()}});
        gelex::bayes::BayesPrior prior{
            gelex::bayes::RandomPrior{make_variance(0.3)},
            std::move(priors),
            gelex::bayes::ResidualPrior{make_variance(0.4)}};
        gelex::BayesState state(model, prior);

        std::mt19937_64 rng{123};
        auto chain = gelex::mcmc::Chain::make(model, prior, state, rng);
        chain.step();

        const auto& block = std::get<gelex::bayes::SingleGeneticBlockState>(
            state.genetics()[0]);
        REQUIRE(block.state().coeffs.allFinite());
        REQUIRE(std::isfinite(block.state().heritability));
    }

    SECTION("per-marker spike-slab")
    {
        std::vector<gelex::bayes::GeneticDesign> genetics;
        genetics.push_back(make_design());
        gelex::BayesModel model{
            Eigen::VectorXd{{1.0, -0.5, 0.25}},
            gelex::FixedDesign::make(3),
            {},
            std::move(genetics)};
        std::vector<gelex::bayes::GeneticPrior> priors;
        priors.emplace_back(
            gelex::bayes::SingleGeneticPrior{
                gelex::bayes::SinglePerMarkerSpikeSlabGaussianPrior{
                    gelex::GeneticMode::A,
                    gelex::bayes::PerMarkerVariance{make_variance(0.1)},
                    make_proportion()}});
        gelex::bayes::BayesPrior prior{
            gelex::bayes::RandomPrior{make_variance(0.3)},
            std::move(priors),
            gelex::bayes::ResidualPrior{make_variance(0.4)}};
        gelex::BayesState state(model, prior);

        std::mt19937_64 rng{123};
        auto chain = gelex::mcmc::Chain::make(model, prior, state, rng);
        chain.step();

        const auto& block = std::get<gelex::bayes::SingleGeneticBlockState>(
            state.genetics()[0]);
        REQUIRE(block.state().coeffs.allFinite());
        REQUIRE(std::isfinite(block.state().heritability));
    }

    SECTION("scaled mixture")
    {
        std::vector<gelex::bayes::GeneticDesign> genetics;
        genetics.push_back(make_design());
        gelex::BayesModel model{
            Eigen::VectorXd{{1.0, -0.5, 0.25}},
            gelex::FixedDesign::make(3),
            {},
            std::move(genetics)};
        std::vector<gelex::bayes::GeneticPrior> priors;
        priors.emplace_back(
            gelex::bayes::SingleGeneticPrior{
                gelex::bayes::SingleScaledMixtureGaussianPrior{
                    gelex::GeneticMode::A,
                    gelex::bayes::SharedMarkerVariance{make_variance(0.1)},
                    Eigen::VectorXd{{0.0, 2.0}},
                    make_proportion()}});
        gelex::bayes::BayesPrior prior{
            gelex::bayes::RandomPrior{make_variance(0.3)},
            std::move(priors),
            gelex::bayes::ResidualPrior{make_variance(0.4)}};
        gelex::BayesState state(model, prior);

        std::mt19937_64 rng{123};
        auto chain = gelex::mcmc::Chain::make(model, prior, state, rng);
        chain.step();

        const auto& block = std::get<gelex::bayes::SingleGeneticBlockState>(
            state.genetics()[0]);
        REQUIRE(block.state().coeffs.allFinite());
        REQUIRE(std::isfinite(block.state().heritability));
    }
}

TEST_CASE(
    "Chain::make maps joint genetic priors to MCMC steps",
    "[mcmc][chain]")
{
    std::vector<gelex::bayes::GeneticDesign> genetics;
    genetics.push_back(make_design(gelex::GeneticMode::A));
    genetics.push_back(make_design(gelex::GeneticMode::D));
    gelex::BayesModel model{
        Eigen::VectorXd{{1.0, -0.5, 0.25}},
        gelex::FixedDesign::make(3),
        {},
        std::move(genetics)};
    std::vector<gelex::bayes::GeneticPrior> priors;
    priors.emplace_back(
        gelex::bayes::JointGeneticPrior{gelex::bayes::JointGaussianMixturePrior{
            gelex::bayes::JointSharedMarkerVariance{std::array{
                gelex::bayes::SharedMarkerVariance{make_variance(0.1)},
                gelex::bayes::SharedMarkerVariance{make_variance(0.2)}}},
            make_proportion(4)}});
    gelex::bayes::BayesPrior prior{
        gelex::bayes::RandomPrior{make_variance(0.3)},
        std::move(priors),
        gelex::bayes::ResidualPrior{make_variance(0.4)}};
    gelex::BayesState state(model, prior);

    std::mt19937_64 rng{123};
    auto chain = gelex::mcmc::Chain::make(model, prior, state, rng);
    chain.step();

    const auto& block
        = std::get<gelex::bayes::JointGeneticBlockState>(state.genetics()[0]);
    REQUIRE(block.state(gelex::GeneticMode::A).coeffs.allFinite());
    REQUIRE(block.state(gelex::GeneticMode::D).coeffs.allFinite());
    REQUIRE(std::isfinite(block.state(gelex::GeneticMode::A).heritability));
    REQUIRE(std::isfinite(block.state(gelex::GeneticMode::D).heritability));
}

TEST_CASE("Chain::make maps half normal mixture step", "[mcmc][chain]")
{
    std::vector<gelex::bayes::GeneticDesign> genetics;
    genetics.push_back(make_design(gelex::GeneticMode::A));
    genetics.push_back(make_design(gelex::GeneticMode::D));
    gelex::BayesModel model{
        Eigen::VectorXd{{1.0, -0.5, 0.25}},
        gelex::FixedDesign::make(3),
        {},
        std::move(genetics)};
    std::vector<gelex::bayes::GeneticPrior> priors;
    priors.emplace_back(
        gelex::bayes::JointGeneticPrior{
            gelex::bayes::JointHalfNormalMixturePrior{
                gelex::bayes::JointSharedMarkerVariance{std::array{
                    gelex::bayes::SharedMarkerVariance{make_variance(0.1)},
                    gelex::bayes::SharedMarkerVariance{make_variance(0.2)}}},
                make_proportion(4),
                gelex::bayes::ProbabilityParameter{
                    0.5, gelex::bayes::BetaPrior{1.0, 1.0}}}});
    gelex::bayes::BayesPrior prior{
        gelex::bayes::RandomPrior{make_variance(0.3)},
        std::move(priors),
        gelex::bayes::ResidualPrior{make_variance(0.4)}};
    gelex::BayesState state(model, prior);

    std::mt19937_64 rng{123};
    auto chain = gelex::mcmc::Chain::make(model, prior, state, rng);
    chain.step();

    const auto& block
        = std::get<gelex::bayes::JointGeneticBlockState>(state.genetics()[0]);
    REQUIRE(block.state(gelex::GeneticMode::A).coeffs.allFinite());
    REQUIRE(block.state(gelex::GeneticMode::D).coeffs.allFinite());
    REQUIRE(std::isfinite(block.state(gelex::GeneticMode::A).heritability));
    REQUIRE(std::isfinite(block.state(gelex::GeneticMode::D).heritability));
}

TEST_CASE("Solver::run collects single genetic samples", "[mcmc][solver]")
{
    std::vector<gelex::bayes::GeneticDesign> genetics;
    genetics.push_back(make_design());
    gelex::BayesModel model{
        Eigen::VectorXd{{1.0, -0.5, 0.25}},
        gelex::FixedDesign::make(3),
        {},
        std::move(genetics)};
    std::vector<gelex::bayes::GeneticPrior> priors;
    priors.emplace_back(
        gelex::bayes::SingleGeneticPrior{
            gelex::bayes::SingleSharedGaussianPrior{
                gelex::GeneticMode::A,
                gelex::bayes::SharedMarkerVariance{make_variance(0.1)}}});
    gelex::bayes::BayesPrior prior{
        gelex::bayes::RandomPrior{make_variance(0.3)},
        std::move(priors),
        gelex::bayes::ResidualPrior{make_variance(0.4)}};

    gelex::mcmc::Params params{
        .n_iters = 4, .n_burn_in = 1, .n_thin = 1, .checkpoint_step = 0};
    gelex::mcmc::Solver solver{params};

    int progress_count = 0;
    bool done = false;
    gelex::MCMCObserver observer = [&](const gelex::MCMCEvent& event)
    {
        if (const auto* progress
            = std::get_if<gelex::MCMCProgressEvent>(&event))
        {
            if (progress->done)
            {
                done = true;
            }
            else
            {
                ++progress_count;
                REQUIRE(progress->current <= progress->total);
                REQUIRE(progress->state != nullptr);
            }
        }
    };

    auto result = solver.run(model, std::move(prior), 123, observer);

    REQUIRE(progress_count == params.n_iters);
    REQUIRE(done);
    REQUIRE(result.samples_collected() == params.n_records());

    const auto& fixed = std::get<gelex::stats::RunningStatsResult>(
        result.get("state/fixed/coeffs"));
    REQUIRE(fixed.mean.size() == 1);
    REQUIRE(fixed.mean.allFinite());

    const auto& residual = std::get<gelex::stats::RunningStatsResult>(
        result.get("state/residual/variance"));
    REQUIRE(residual.mean.size() == 1);
    REQUIRE(residual.mean.allFinite());

    const auto& coeffs = std::get<gelex::stats::RunningStatsResult>(
        result.get("state/genetic_0/single/A/genetic/coeffs"));
    REQUIRE(coeffs.mean.size() == 2);
    REQUIRE(coeffs.mean.allFinite());

    const auto& pve = std::get<gelex::stats::RunningStatsResult>(
        result.get("state/genetic_0/single/A/genetic/pve"));
    REQUIRE(pve.mean.size() == 2);
    REQUIRE(pve.mean.allFinite());

    const auto& variance = std::get<gelex::stats::RunningStatsResult>(
        result.get("state/genetic_0/single/A/genetic/variance"));
    REQUIRE(variance.mean.size() == 1);
    REQUIRE(variance.mean.allFinite());

    const auto& heritability = std::get<gelex::stats::RunningStatsResult>(
        result.get("state/genetic_0/single/A/genetic/heritability"));
    REQUIRE(heritability.mean.size() == 1);
    REQUIRE(heritability.mean.allFinite());
}

TEST_CASE("Solver::run collects joint genetic samples", "[mcmc][solver]")
{
    std::vector<gelex::bayes::GeneticDesign> genetics;
    genetics.push_back(make_design(gelex::GeneticMode::A));
    genetics.push_back(make_design(gelex::GeneticMode::D));
    gelex::BayesModel model{
        Eigen::VectorXd{{1.0, -0.5, 0.25}},
        gelex::FixedDesign::make(3),
        {},
        std::move(genetics)};
    std::vector<gelex::bayes::GeneticPrior> priors;
    priors.emplace_back(
        gelex::bayes::JointGeneticPrior{gelex::bayes::JointGaussianMixturePrior{
            gelex::bayes::JointSharedMarkerVariance{std::array{
                gelex::bayes::SharedMarkerVariance{make_variance(0.1)},
                gelex::bayes::SharedMarkerVariance{make_variance(0.2)}}},
            make_proportion(4)}});
    gelex::bayes::BayesPrior prior{
        gelex::bayes::RandomPrior{make_variance(0.3)},
        std::move(priors),
        gelex::bayes::ResidualPrior{make_variance(0.4)}};

    gelex::mcmc::Params params{
        .n_iters = 4, .n_burn_in = 1, .n_thin = 1, .checkpoint_step = 0};
    gelex::mcmc::Solver solver{params};
    auto result = solver.run(model, std::move(prior), 123);

    const auto& additive_coeffs = std::get<gelex::stats::RunningStatsResult>(
        result.get("state/genetic_0/joint/A/genetic/coeffs"));
    const auto& dominance_coeffs = std::get<gelex::stats::RunningStatsResult>(
        result.get("state/genetic_0/joint/D/genetic/coeffs"));
    REQUIRE(additive_coeffs.mean.size() == 2);
    REQUIRE(dominance_coeffs.mean.size() == 2);
    REQUIRE(additive_coeffs.mean.allFinite());
    REQUIRE(dominance_coeffs.mean.allFinite());

    const auto& additive_pve = std::get<gelex::stats::RunningStatsResult>(
        result.get("state/genetic_0/joint/A/genetic/pve"));
    const auto& dominance_pve = std::get<gelex::stats::RunningStatsResult>(
        result.get("state/genetic_0/joint/D/genetic/pve"));
    REQUIRE(additive_pve.mean.size() == 2);
    REQUIRE(dominance_pve.mean.size() == 2);
    REQUIRE(additive_pve.mean.allFinite());
    REQUIRE(dominance_pve.mean.allFinite());

    const auto& additive_pip
        = std::get<gelex::stats::RunningStatsResult>(result.get(
            "state/genetic_0/joint/prior_state/"
            "joint_mixture_gaussian/mixture/A/pip"));
    const auto& dominance_pip
        = std::get<gelex::stats::RunningStatsResult>(result.get(
            "state/genetic_0/joint/prior_state/"
            "joint_mixture_gaussian/mixture/D/pip"));
    REQUIRE(additive_pip.mean.size() == 2);
    REQUIRE(dominance_pip.mean.size() == 2);
    REQUIRE(additive_pip.mean.allFinite());
    REQUIRE(dominance_pip.mean.allFinite());

    const auto& additive_variance = std::get<gelex::stats::RunningStatsResult>(
        result.get("state/genetic_0/joint/A/genetic/variance"));
    const auto& dominance_variance = std::get<gelex::stats::RunningStatsResult>(
        result.get("state/genetic_0/joint/D/genetic/variance"));
    REQUIRE(additive_variance.mean.allFinite());
    REQUIRE(dominance_variance.mean.allFinite());

    const auto& additive_heritability
        = std::get<gelex::stats::RunningStatsResult>(
            result.get("state/genetic_0/joint/A/genetic/heritability"));
    const auto& dominance_heritability
        = std::get<gelex::stats::RunningStatsResult>(
            result.get("state/genetic_0/joint/D/genetic/heritability"));
    REQUIRE(additive_heritability.mean.allFinite());
    REQUIRE(dominance_heritability.mean.allFinite());
}

TEST_CASE("MCMC command dataflow writes solver outputs", "[mcmc][solver]")
{
    std::vector<gelex::bayes::GeneticDesign> genetics;
    genetics.push_back(make_design());
    gelex::BayesModel model{
        Eigen::VectorXd{{1.0, -0.5, 0.25}},
        gelex::FixedDesign::make(3),
        {},
        std::move(genetics)};
    std::vector<gelex::bayes::GeneticPrior> priors;
    priors.emplace_back(
        gelex::bayes::SingleGeneticPrior{
            gelex::bayes::SingleSharedGaussianPrior{
                gelex::GeneticMode::A,
                gelex::bayes::SharedMarkerVariance{make_variance(0.1)}}});
    gelex::bayes::BayesPrior prior{
        gelex::bayes::RandomPrior{make_variance(0.3)},
        std::move(priors),
        gelex::bayes::ResidualPrior{make_variance(0.4)}};

    gelex::mcmc::Params params{
        .n_iters = 4, .n_burn_in = 1, .n_thin = 1, .checkpoint_step = 2};
    gelex::test::FileFixture files;
    const auto prefix = files.get_test_dir() / "mcmc_solver_output_test";
    const auto bfile_prefix = files.get_test_dir() / "mcmc_solver_output_geno";
    (void)files.create_named_text_file(
        "mcmc_solver_output_geno.bim",
        "1\trs1\t0\t100\tA\tG\n"
        "1\trs2\t0\t200\tC\tT\n");
    gelex::mcmc::Solver solver{
        params, prefix.string() + ".draws", prefix.string()};

    bool done = false;
    int checkpoint_saved = 0;
    gelex::MCMCObserver observer = [&](const gelex::MCMCEvent& event)
    {
        if (const auto* progress
            = std::get_if<gelex::MCMCProgressEvent>(&event))
        {
            if (progress->done)
            {
                done = true;
            }
            else
            {
                REQUIRE(progress->state != nullptr);
                REQUIRE(
                    prior.genetics().size()
                    == progress->state->genetics().size());
            }
        }
        if (std::get_if<gelex::MCMCCheckpointSavedEvent>(&event) != nullptr)
        {
            ++checkpoint_saved;
        }
    };

    auto result = solver.run(model, prior, 123, observer);
    gelex::mcmc::write_params(result, prefix.string());
    gelex::mcmc::write_summary(result, prefix.string());
    gelex::mcmc::write_snp_eff(
        result, model, bfile_prefix.string() + ".bim", prefix.string());

    REQUIRE(done);
    REQUIRE(checkpoint_saved == 2);
    REQUIRE(result.samples_collected() == params.n_records());
    auto summary_path = prefix;
    summary_path += ".summary";
    REQUIRE(std::filesystem::exists(summary_path));
    auto params_path = prefix;
    params_path += ".param";
    REQUIRE(std::filesystem::exists(params_path));
    auto snp_path = prefix;
    snp_path += ".snpeff";
    REQUIRE(std::filesystem::exists(snp_path));
    auto draws_path = prefix;
    draws_path += ".draws";
    REQUIRE(std::filesystem::exists(draws_path));
    auto checkpoint_path = prefix;
    checkpoint_path += ".ckpt";
    REQUIRE(std::filesystem::exists(checkpoint_path));
    gelex::io::BinaryReader draws(draws_path.string());
    const auto fixed_draws = draws.to_map<double>("state/fixed/coeffs");
    REQUIRE(fixed_draws.rows() == 1);
    REQUIRE(fixed_draws.cols() == params.n_records());
}

TEST_CASE("Solver rejects invalid MCMC params", "[mcmc][solver]")
{
    SECTION("n_iters must be positive")
    {
        REQUIRE_THROWS_AS(
            gelex::mcmc::Solver(
                gelex::mcmc::Params{
                    .n_iters = 0,
                    .n_burn_in = 0,
                    .n_thin = 1,
                    .checkpoint_step = 0}),
            gelex::GelexException);
    }

    SECTION("burn-in must be in range")
    {
        REQUIRE_THROWS_AS(
            gelex::mcmc::Solver(
                gelex::mcmc::Params{
                    .n_iters = 4,
                    .n_burn_in = 4,
                    .n_thin = 1,
                    .checkpoint_step = 0}),
            gelex::GelexException);
    }

    SECTION("thin must be positive")
    {
        REQUIRE_THROWS_AS(
            gelex::mcmc::Solver(
                gelex::mcmc::Params{
                    .n_iters = 4,
                    .n_burn_in = 1,
                    .n_thin = 0,
                    .checkpoint_step = 0}),
            gelex::GelexException);
    }

    SECTION("thin must divide retained iterations")
    {
        REQUIRE_THROWS_AS(
            gelex::mcmc::Solver(
                gelex::mcmc::Params{
                    .n_iters = 5,
                    .n_burn_in = 1,
                    .n_thin = 3,
                    .checkpoint_step = 0}),
            gelex::GelexException);
    }

    SECTION("checkpoint step must be non-negative")
    {
        REQUIRE_THROWS_AS(
            gelex::mcmc::Solver(
                gelex::mcmc::Params{
                    .n_iters = 4,
                    .n_burn_in = 1,
                    .n_thin = 1,
                    .checkpoint_step = -1}),
            gelex::GelexException);
    }
}

TEST_CASE(
    "Solver::run_from starts from checkpointed state",
    "[mcmc][solver][checkpoint]")
{
    std::vector<gelex::bayes::GeneticDesign> genetics;
    genetics.push_back(make_design());
    gelex::BayesModel model{
        Eigen::VectorXd{{1.0, -0.5, 0.25}},
        gelex::FixedDesign::make(3),
        {},
        std::move(genetics)};
    gelex::test::FileFixture files;
    const auto continuous_prefix
        = (files.get_test_dir() / "solver_continuous").string();
    const auto first_prefix = (files.get_test_dir() / "solver_first").string();
    const auto from_prefix = (files.get_test_dir() / "solver_from").string();

    std::vector<gelex::bayes::GeneticPrior> continuous_priors;
    continuous_priors.emplace_back(
        gelex::bayes::SingleGeneticPrior{
            gelex::bayes::SingleSharedGaussianPrior{
                gelex::GeneticMode::A,
                gelex::bayes::SharedMarkerVariance{make_variance(0.1)}}});
    gelex::bayes::BayesPrior continuous_prior{
        gelex::bayes::RandomPrior{make_variance(0.3)},
        std::move(continuous_priors),
        gelex::bayes::ResidualPrior{make_variance(0.4)}};
    gelex::mcmc::Params continuous_params{
        .n_iters = 8, .n_burn_in = 0, .n_thin = 1, .checkpoint_step = 8};
    gelex::mcmc::Solver continuous_solver{
        continuous_params, "", continuous_prefix};
    continuous_solver.run(model, continuous_prior, 123);

    std::vector<gelex::bayes::GeneticPrior> first_priors;
    first_priors.emplace_back(
        gelex::bayes::SingleGeneticPrior{
            gelex::bayes::SingleSharedGaussianPrior{
                gelex::GeneticMode::A,
                gelex::bayes::SharedMarkerVariance{make_variance(0.1)}}});
    gelex::bayes::BayesPrior first_prior{
        gelex::bayes::RandomPrior{make_variance(0.3)},
        std::move(first_priors),
        gelex::bayes::ResidualPrior{make_variance(0.4)}};
    gelex::mcmc::Params first_params{
        .n_iters = 4, .n_burn_in = 0, .n_thin = 1, .checkpoint_step = 4};
    gelex::mcmc::Solver first_solver{first_params, "", first_prefix};
    first_solver.run(model, first_prior, 123);

    std::vector<gelex::bayes::GeneticPrior> from_priors;
    from_priors.emplace_back(
        gelex::bayes::SingleGeneticPrior{
            gelex::bayes::SingleSharedGaussianPrior{
                gelex::GeneticMode::A,
                gelex::bayes::SharedMarkerVariance{make_variance(0.1)}}});
    gelex::bayes::BayesPrior from_prior{
        gelex::bayes::RandomPrior{make_variance(0.3)},
        std::move(from_priors),
        gelex::bayes::ResidualPrior{make_variance(0.4)}};
    gelex::mcmc::Params from_params{
        .n_iters = 4, .n_burn_in = 0, .n_thin = 1, .checkpoint_step = 4};
    gelex::mcmc::Solver from_solver{from_params, "", from_prefix};
    const auto from_result = from_solver.run_from(
        model, from_prior, first_prefix + ".ckpt");

    REQUIRE(from_result.samples_collected() == from_params.n_records());

    gelex::io::BinaryReader continuous{continuous_prefix + ".ckpt"};
    gelex::io::BinaryReader from_checkpoint{from_prefix + ".ckpt"};
    CHECK(continuous.to_map<double>("state/fixed/coeffs")
              .isApprox(
                  from_checkpoint.to_map<double>("state/fixed/coeffs"), 0.0));
    CHECK(continuous.to_map<double>("state/genetic_0/single/A/genetic/coeffs")
              .isApprox(
                  from_checkpoint.to_map<double>(
                      "state/genetic_0/single/A/genetic/coeffs"),
                  0.0));
    CHECK(continuous.to_map<double>("state/genetic_0/single/A/genetic/u")
              .isApprox(
                  from_checkpoint.to_map<double>(
                      "state/genetic_0/single/A/genetic/u"),
                  0.0));
    CHECK(continuous.to_map<double>("state/genetic_0/single/A/genetic/variance")
              .isApprox(
                  from_checkpoint.to_map<double>(
                      "state/genetic_0/single/A/genetic/variance"),
                  0.0));
    CHECK(continuous
              .to_map<double>("state/genetic_0/single/A/prior_state/"
                              "shared_gaussian/variance")
              .isApprox(
                  from_checkpoint.to_map<double>(
                      "state/genetic_0/single/A/prior_state/"
                      "shared_gaussian/variance"),
                  0.0));
    CHECK(continuous.to_map<double>("state/residual/y_adj")
              .isApprox(
                  from_checkpoint.to_map<double>("state/residual/y_adj"), 0.0));
    CHECK(continuous.to_map<double>("state/residual/variance")
              .isApprox(
                  from_checkpoint.to_map<double>("state/residual/variance"),
                  0.0));

    const auto continuous_rng = continuous.to_strings("rng_state");
    const auto from_rng = from_checkpoint.to_strings("rng_state");
    REQUIRE(continuous_rng.size() == 1);
    REQUIRE(from_rng.size() == 1);
    CHECK(continuous_rng.front() == from_rng.front());
}
