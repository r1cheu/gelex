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
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <random>
#include <utility>
#include <variant>
#include <vector>

#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/gaussian_prior_state.h"
#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/legacy_prior.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/state.h"
#include "gelex/exception.h"
#include "gelex/io/binary_writer.h"
#include "gelex/io/mcmc_checkpoint.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_mode.h"

#include "compact_genotype_fixture.h"
#include "file_fixture.h"

namespace
{

auto make_variance(double value) -> gelex::bayes::VarianceParameter
{
    return gelex::bayes::VarianceParameter{
        value, gelex::bayes::ScaledInvChiSqPrior{4.0, 1.0}};
}

auto make_model() -> gelex::BayesModel
{
    return gelex::test::make_compact_model(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}},
        Eigen::VectorXd{{1.0, 2.0, 3.0}});
}

auto make_prior() -> gelex::bayes::BayesPrior
{
    std::vector<gelex::bayes::GeneticPrior> genetics;
    genetics.emplace_back(
        gelex::bayes::SingleGeneticPrior{
            gelex::bayes::SingleSharedGaussianPrior{
                gelex::GeneticMode::A,
                gelex::bayes::SharedMarkerVariance{make_variance(0.1)}}});

    return gelex::bayes::BayesPrior{
        gelex::bayes::RandomPrior{make_variance(0.3)},
        std::move(genetics),
        gelex::bayes::ResidualPrior{make_variance(0.4)}};
}

auto make_state() -> gelex::BayesState
{
    auto model = make_model();
    auto prior = make_prior();
    return gelex::BayesState{model, prior};
}

}  // namespace

TEST_CASE(
    "MCMC checkpoint round-trip preserves BayesState fields",
    "[checkpoint]")
{
    gelex::test::FileFixture files;
    const auto prefix = (files.get_test_dir() / "mcmc_state").string();
    const auto checkpoint_path = prefix + ".ckpt";

    auto state = make_state();
    state.fixed().coeffs = Eigen::VectorXd{{1.25}};
    auto& block
        = std::get<gelex::bayes::SingleGeneticBlockState>(state.genetics()[0]);
    block.state().coeffs = Eigen::VectorXd{{0.1, 0.2}};
    block.state().u = Eigen::VectorXd{{-0.5, 0.0, 0.5}};
    block.state().variance = 0.75;
    block.state().heritability = 0.3;
    auto& prior_state = std::get<gelex::bayes::SingleSharedGaussianState>(
        block.prior_state());
    prior_state.variance() = 0.45;
    state.residual().y_adj = Eigen::VectorXd{{0.25, 0.5, 0.75}};
    state.residual().variance = 1.5;

    std::mt19937_64 rng{123};
    rng();
    rng();

    gelex::write_checkpoint(state, rng, prefix);

    REQUIRE(std::filesystem::exists(checkpoint_path));

    auto restored = make_state();
    auto restored_rng = gelex::read_checkpoint(checkpoint_path, restored);
    const auto& restored_block
        = std::get<gelex::bayes::SingleGeneticBlockState>(
            restored.genetics()[0]);
    const auto& restored_prior_state
        = std::get<gelex::bayes::SingleSharedGaussianState>(
            restored_block.prior_state());

    CHECK(restored.fixed().coeffs.isApprox(state.fixed().coeffs));
    CHECK(restored_block.state().coeffs.isApprox(block.state().coeffs));
    CHECK(restored_block.state().u.isApprox(block.state().u));
    CHECK(restored_block.state().variance == block.state().variance);
    CHECK(restored_block.state().heritability == block.state().heritability);
    CHECK(restored_prior_state.variance() == prior_state.variance());
    CHECK(restored.residual().y_adj.isApprox(state.residual().y_adj));
    CHECK(restored.residual().variance == state.residual().variance);
    CHECK(restored_rng() == rng());
}

TEST_CASE("MCMC checkpoint reader rejects field shape mismatch", "[checkpoint]")
{
    gelex::test::FileFixture files;
    const auto checkpoint_path = files.generate_random_file_path(".ckpt");
    gelex::BinaryWriter writer(checkpoint_path.string());
    const Eigen::MatrixXd fixed_coeffs{{1.0}, {2.0}};
    writer.write("state/fixed/coeffs", fixed_coeffs);
    writer.close();

    auto state = make_state();

    REQUIRE_THROWS_AS(
        gelex::read_checkpoint(checkpoint_path, state), gelex::GelexException);
}
