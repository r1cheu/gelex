// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <type_traits>
#include <vector>

#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/mcmc_runner.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/spec.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_reader.h"

#include "compact_genotype_fixture.h"
#include "file_fixture.h"

namespace
{

constexpr auto mode_a = gelex::GeneticModeSet{gelex::GeneticMode::A};

}  // namespace

TEST_CASE(
    "MCMC runner executes typed kernels and reports every iteration",
    "[bayes][mcmc][runner]")
{
    const auto model = gelex::test::make_compact_model(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0}, {0.0, 1.0}},
        Eigen::VectorXd{{1.0, -0.5, 0.25, 2.0}});
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<
            mode_a,
            gelex::HomogeneousModeValues<
                mode_a,
                gelex::GaussianSpec<gelex::VarianceLayout::Pooled>>>::
            defaults(),
        model);
    constexpr int iterations = 4;
    gelex::MCMCRunner runner{iterations, 0, 1};
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "progress.draws";
    std::vector<std::size_t> completed_iterations;
    const auto observer
        = [&](std::size_t current) { completed_iterations.push_back(current); };

    static_assert(
        std::is_void_v<decltype(runner.run(model, prior, path.string()))>);
    runner.run(model, prior, path.string(), 123, observer);

    REQUIRE(completed_iterations == std::vector<std::size_t>{1, 2, 3, 4});
    const gelex::BinaryReader reader{path.string()};
    REQUIRE(reader.to_map<double>("residual/variance").cols() == iterations);
}

TEST_CASE(
    "MCMC runner retains draws after burn-in at the thinning interval",
    "[bayes][mcmc][runner]")
{
    const auto model = gelex::test::make_compact_model(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0}, {0.0, 1.0}},
        Eigen::VectorXd{{1.0, -0.5, 0.25, 2.0}});
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<
            mode_a,
            gelex::HomogeneousModeValues<
                mode_a,
                gelex::GaussianSpec<gelex::VarianceLayout::Pooled>>>::
            defaults(),
        model);
    gelex::test::FileFixture fixture;
    const auto full_path = fixture.get_test_dir() / "full.draws";
    const auto retained_path = fixture.get_test_dir() / "retained.draws";

    gelex::MCMCRunner full_runner{5, 0, 1};
    full_runner.run(model, prior, full_path.string(), 123);
    gelex::MCMCRunner retained_runner{5, 1, 2};
    REQUIRE(retained_runner.draw_count() == 2);
    retained_runner.run(model, prior, retained_path.string(), 123);

    const gelex::BinaryReader full_reader(full_path.string());
    const gelex::BinaryReader retained_reader(retained_path.string());
    const auto full = full_reader.to_map<double>("residual/variance");
    const auto retained = retained_reader.to_map<double>("residual/variance");
    const Eigen::MatrixXd expected{{full(0, 2), full(0, 4)}};
    REQUIRE(retained.isApprox(expected));
}

TEST_CASE("MCMC runner rejects invalid schedules", "[bayes][mcmc][runner]")
{
    REQUIRE_THROWS_AS(gelex::MCMCRunner(0, 0, 1), gelex::GelexException);
    REQUIRE_THROWS_AS(gelex::MCMCRunner(-1, 0, 1), gelex::GelexException);
    REQUIRE_THROWS_AS(gelex::MCMCRunner(4, -1, 1), gelex::GelexException);
    REQUIRE_THROWS_AS(gelex::MCMCRunner(4, 4, 1), gelex::GelexException);
    REQUIRE_THROWS_AS(gelex::MCMCRunner(4, 0, 0), gelex::GelexException);
    REQUIRE_THROWS_AS(gelex::MCMCRunner(4, 0, -1), gelex::GelexException);
    REQUIRE_THROWS_AS(gelex::MCMCRunner(4, 1, 2), gelex::GelexException);
}
