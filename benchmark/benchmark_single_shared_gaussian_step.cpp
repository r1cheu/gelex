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

#include <filesystem>
#include <random>
#include <string_view>
#include <utility>
#include <variant>

#include <Eigen/Core>

#include <nanobench.h>

#include <catch2/catch_test_macros.hpp>

#include "gelex/algo/mcmc/steps/single_genetic_step.h"
#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/gaussian_prior_state.h"
#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/parameter/distributions.h"
#include "gelex/bayes/parameter/values.h"
#include "gelex/bayes/state.h"
#include "gelex/data/genotype_reader.h"
#include "gelex/types/genetic_mode.h"

namespace
{

constexpr std::string_view kGbinPath
    = "/home/rlchen/tmp/gelex_single_shared_gaussian_step.gbin";

}  // namespace

TEST_CASE(
    "SingleSharedGaussianStep throughput",
    "[!benchmark][mcmc][genetic][single_shared_gaussian]")
{
    REQUIRE(std::filesystem::exists(kGbinPath));

    auto genotype = gelex::genotype::GenotypeReader::read(
        std::filesystem::path{kGbinPath}, gelex::GeneticMode::A);
    gelex::bayes::GeneticDesign design{
        gelex::GeneticMode::A, std::move(genotype)};
    gelex::bayes::SingleGeneticPrior prior{
        gelex::bayes::SingleSharedGaussianPrior{
            gelex::GeneticMode::A,
            gelex::bayes::SharedMarkerVariance{gelex::bayes::VarianceParameter{
                0.1, gelex::bayes::ScaledInvChiSqPrior{4.0, 1.0}}}}};
    gelex::bayes::SingleGeneticBlockState block{design, prior};
    gelex::bayes::ResidualState residual{
        .y_adj = Eigen::VectorXd::Zero(design.X.rows()),
        .variance = 1.0,
        .old_y_adj = Eigen::VectorXd::Zero(design.X.rows())};
    for (Eigen::Index i = 0; i < residual.y_adj.size(); ++i)
    {
        residual.y_adj(i) = 0.01 * static_cast<double>((i % 251) - 125);
    }
    auto& shared_state = std::get<gelex::bayes::SingleSharedGaussianState>(
        block.prior_state());

    std::mt19937_64 rng{42};
    gelex::mcmc::SingleSharedGaussianStep step{
        design, prior, block, residual, rng};

    ankerl::nanobench::Bench b;
    b.title("SingleSharedGaussianStep")
        .unit("marker")
        .batch(static_cast<double>(design.X.cols()))
        .warmup(3)
        .minEpochIterations(60);

    b.run(
        "step",
        [&]
        {
            step.step();
            ankerl::nanobench::doNotOptimizeAway(block.state().coeffs.data());
            ankerl::nanobench::doNotOptimizeAway(residual.y_adj.data());
            ankerl::nanobench::doNotOptimizeAway(shared_state.variance());
        });
}
