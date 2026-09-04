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
#include <cstddef>
#include <type_traits>
#include <vector>

#include "gelex/bayes/draws.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/mcmc_runner.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_reader.h"

#include "compact_genotype_fixture.h"
#include "file_fixture.h"

namespace
{

constexpr auto mode_a = gelex::GeneticModeSet{gelex::GeneticMode::A};
using Family = gelex::GaussianFamily<gelex::VarianceLayout::Pooled>;

}  // namespace

TEST_CASE(
    "MCMC runner executes typed kernels and reports every iteration",
    "[bayes][mcmc][runner]")
{
    const auto model = gelex::test::make_compact_model(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0}, {0.0, 1.0}},
        Eigen::VectorXd{{1.0, -0.5, 0.25, 2.0}});
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_a, Family>::defaults(), model);
    constexpr int iterations = 4;
    gelex::MCMCRunner runner{iterations, 0, 1};
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "progress.draws";
    auto draws
        = gelex::make_draws(prior, model, path.string(), runner.draw_count());
    std::vector<std::size_t> completed_iterations;
    const auto observer
        = [&](std::size_t current) { completed_iterations.push_back(current); };

    static_assert(std::is_void_v<decltype(runner.run(model, prior, draws))>);
    runner.run(model, prior, draws, 123, observer);

    REQUIRE(completed_iterations == std::vector<std::size_t>{1, 2, 3, 4});
}

TEST_CASE(
    "MCMC runner retains draws after burn-in at the thinning interval",
    "[bayes][mcmc][runner]")
{
    const auto model = gelex::test::make_compact_model(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0}, {0.0, 1.0}},
        Eigen::VectorXd{{1.0, -0.5, 0.25, 2.0}});
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_a, Family>::defaults(), model);
    gelex::test::FileFixture fixture;
    const auto full_path = fixture.get_test_dir() / "full.draws";
    const auto retained_path = fixture.get_test_dir() / "retained.draws";

    {
        gelex::MCMCRunner runner{5, 0, 1};
        auto draws = gelex::make_draws(
            prior, model, full_path.string(), runner.draw_count());
        runner.run(model, prior, draws, 123);
    }
    {
        gelex::MCMCRunner runner{5, 1, 2};
        REQUIRE(runner.draw_count() == 2);
        auto draws = gelex::make_draws(
            prior, model, retained_path.string(), runner.draw_count());
        runner.run(model, prior, draws, 123);
    }

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
