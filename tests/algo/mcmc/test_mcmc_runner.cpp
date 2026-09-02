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

#include "gelex/algo/mcmc/mcmc_runner.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/semantic_method.h"
#include "gelex/exception.h"
#include "gelex/types/genetic_mode.h"

#include "../../compact_genotype_fixture.h"

namespace
{

constexpr auto mode_a = gelex::GeneticModeSet{gelex::GeneticMode::A};
using Method = gelex::GaussianMethod<gelex::VarianceLayout::Pooled>;

}  // namespace

TEST_CASE(
    "MCMC runner executes typed kernels and reports every iteration",
    "[algo][mcmc][runner]")
{
    const auto model = gelex::test::make_compact_model(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0}, {0.0, 1.0}},
        Eigen::VectorXd{{1.0, -0.5, 0.25, 2.0}});
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_a, Method>::defaults(), model);
    constexpr int iterations = 4;
    gelex::MCMCRunner runner{iterations};
    std::vector<std::size_t> completed_iterations;
    bool done = false;
    gelex::MCMCObserver observer = [&](const gelex::MCMCProgressEvent& progress)
    {
        REQUIRE(progress.total == static_cast<std::size_t>(iterations));
        if (progress.done)
        {
            done = true;
        }
        else
        {
            completed_iterations.push_back(progress.current);
        }
    };

    static_assert(std::is_void_v<decltype(runner.run(model, prior))>);
    runner.run(model, prior, 123, observer);

    REQUIRE(completed_iterations == std::vector<std::size_t>{1, 2, 3, 4});
    REQUIRE(done);
}

TEST_CASE("MCMC runner rejects invalid iterations", "[algo][mcmc][runner]")
{
    REQUIRE_THROWS_AS(gelex::MCMCRunner(0), gelex::GelexException);
    REQUIRE_THROWS_AS(gelex::MCMCRunner(-1), gelex::GelexException);
}
