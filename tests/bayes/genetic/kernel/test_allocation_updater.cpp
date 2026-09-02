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
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>

#include "gelex/bayes/genetic/kernel/allocation_updater.h"
#include "gelex/bayes/genetic/parameter.h"
#include "gelex/bayes/genetic_family.h"

using Catch::Approx;

TEST_CASE(
    "categorical allocation builds and draws the class posterior",
    "[bayes][kernel][allocation]")
{
    constexpr std::array prior_probabilities{0.2, 0.3, 0.5};
    gelex::detail::CategoricalAllocationUpdater<
        gelex::MixtureWeightUpdate::Disabled,
        prior_probabilities.size()>
        allocation{gelex::FixedParameter{prior_probabilities}};
    allocation.begin_sweep(prior_probabilities);

    const auto posterior
        = allocation.posterior({0.0, std::log(2.0), std::log(4.0)});
    REQUIRE(std::exp(posterior.log_marginal_kernel) == Approx(2.8));

    constexpr std::size_t draw_count = 20'000;
    std::array<std::size_t, prior_probabilities.size()> counts{};
    std::mt19937_64 rng{123};
    for (std::size_t draw_index = 0; draw_index < draw_count; ++draw_index)
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
        ++counts[allocation.draw(posterior, rng)];
    }

    REQUIRE(
        static_cast<double>(counts[0]) / draw_count
        == Approx(0.2 / 2.8).margin(0.01));
    REQUIRE(
        static_cast<double>(counts[1]) / draw_count
        == Approx(0.6 / 2.8).margin(0.01));
    REQUIRE(
        static_cast<double>(counts[2]) / draw_count
        == Approx(2.0 / 2.8).margin(0.01));

    auto updated_probabilities = prior_probabilities;
    allocation.update(updated_probabilities, rng);
    REQUIRE(updated_probabilities == prior_probabilities);
}

TEST_CASE(
    "allocation posterior remains stable under extreme log weights",
    "[bayes][kernel][allocation]")
{
    constexpr std::array prior_probabilities{0.25, 0.75};
    gelex::detail::CategoricalAllocationUpdater<
        gelex::MixtureWeightUpdate::Disabled,
        prior_probabilities.size()>
        allocation{gelex::FixedParameter{prior_probabilities}};
    allocation.begin_sweep(prior_probabilities);

    const auto posterior = allocation.posterior({1'000.0, 1'002.0});
    const double expected = 1'000.0 + std::log(0.25 + (0.75 * std::exp(2.0)));
    REQUIRE(posterior.log_marginal_kernel == Approx(expected));
    REQUIRE(std::isfinite(posterior.total_weight));
}

TEST_CASE(
    "binary allocation records sampled outcomes for the weight update",
    "[bayes][kernel][allocation]")
{
    gelex::detail::BinaryAllocationUpdater<gelex::MixtureWeightUpdate::Enabled>
        allocation{gelex::SampledParameter<double, gelex::BetaHyperPrior>{
            .initial = 0.5, .hyperprior = {1.0, 1.0}}};
    allocation.begin_sweep(0.5);
    const auto posterior
        = allocation.posterior({-std::numeric_limits<double>::infinity(), 0.0});

    std::mt19937_64 rng{321};
    for (std::size_t draw_index = 0; draw_index < 100; ++draw_index)
    {
        REQUIRE(allocation.draw(posterior, rng));
    }

    double probability = 0.5;
    allocation.update(probability, rng);
    REQUIRE(probability > 0.9);
    REQUIRE(probability < 1.0);
}
