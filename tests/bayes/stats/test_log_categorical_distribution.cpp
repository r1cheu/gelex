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
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <type_traits>

#include "gelex/bayes/stats/log_categorical_distribution.h"

using Catch::Approx;
using gelex::LogCategoricalDistribution;
using gelex::make_log_weights;
using gelex::make_mixture_posterior_weights;

template <std::size_t K>
concept ValidLogCategoricalDistribution
    = requires { typename LogCategoricalDistribution<K>; };

static_assert(ValidLogCategoricalDistribution<1>);
static_assert(!ValidLogCategoricalDistribution<0>);

TEST_CASE(
    "Log weights are derived from probabilities",
    "[bayes][stats][log_categorical_distribution]")
{
    const auto log_weight_values = make_log_weights(std::array{0.2, 0.3, 0.5});
    const Eigen::Map<const Eigen::Vector3d> log_weights{
        log_weight_values.data()};
    const Eigen::Vector3d expected{
        {std::log(0.2), std::log(0.3), std::log(0.5)}};

    REQUIRE(log_weights.isApprox(expected));
}

TEST_CASE(
    "LogCategoricalDistribution exposes prepared parameters",
    "[bayes][stats][log_categorical_distribution]")
{
    using Distribution = LogCategoricalDistribution<3>;
    using Parameters = Distribution::param_type;
    static_assert(std::is_same_v<Parameters::distribution_type, Distribution>);

    const Parameters default_parameters;
    const auto default_probability_values = default_parameters.probabilities();
    const Eigen::Map<const Eigen::Vector3d> default_probabilities{
        default_probability_values.data()};
    const Eigen::Vector3d expected_default
        = Eigen::Vector3d::Constant(1.0 / 3.0);
    REQUIRE(default_probabilities.isApprox(expected_default));

    const Parameters parameters{Distribution::log_weights_type{
        std::log(1.0), std::log(2.0), std::log(3.0)}};
    const auto probability_values = parameters.probabilities();
    const Eigen::Map<const Eigen::Vector3d> probabilities{
        probability_values.data()};
    const Eigen::Vector3d expected{{1.0 / 6.0, 2.0 / 6.0, 3.0 / 6.0}};
    REQUIRE(probabilities.isApprox(expected));
}

TEST_CASE(
    "Mixture posterior weights combine mixture weights and component integrals",
    "[bayes][stats][log_categorical_distribution]")
{
    const auto parameters = make_mixture_posterior_weights(
        std::array{std::log(0.25), std::log(0.75)},
        std::array{std::log(2.0), std::log(4.0)});
    const auto probability_values = parameters.probabilities();
    const Eigen::Map<const Eigen::Vector2d> probabilities{
        probability_values.data()};
    const Eigen::Vector2d expected{{1.0 / 7.0, 6.0 / 7.0}};

    REQUIRE(probabilities.isApprox(expected));
}

TEST_CASE(
    "LogCategoricalDistribution follows the STL parameter API",
    "[bayes][stats][log_categorical_distribution]")
{
    using Distribution = LogCategoricalDistribution<3>;
    const Distribution::log_weights_type log_weights{0.0, 1.0, 2.0};
    const Distribution::param_type parameters{log_weights};

    Distribution distribution;
    distribution.param(parameters);
    REQUIRE(distribution.param() == parameters);
    REQUIRE(distribution == Distribution{parameters});
    REQUIRE(distribution == Distribution{log_weights});
    REQUIRE(distribution.probabilities() == parameters.probabilities());
    REQUIRE(distribution.min() == 0);
    REQUIRE(distribution.max() == 2);

    distribution.reset();
    REQUIRE(distribution.param() == parameters);

    std::mt19937_64 stored_rng{481};
    std::mt19937_64 supplied_rng{481};
    Distribution default_distribution;
    for (std::size_t draw_index = 0; draw_index < 100; ++draw_index)
    {
        REQUIRE(
            distribution(stored_rng)
            == default_distribution(supplied_rng, parameters));
    }
}

TEST_CASE(
    "LogCategoricalDistribution samples normalized log weights",
    "[bayes][stats][log_categorical_distribution]")
{
    constexpr std::size_t draw_count = 20'000;
    LogCategoricalDistribution<3> distribution{
        std::array{std::log(0.2), std::log(0.3), std::log(0.5)}};
    std::array<std::size_t, 3> counts{};
    std::mt19937_64 rng{123};

    for (std::size_t draw_index = 0; draw_index < draw_count; ++draw_index)
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
        ++counts[distribution(rng)];
    }

    REQUIRE(
        static_cast<double>(counts[0]) / draw_count
        == Approx(0.2).margin(0.01));
    REQUIRE(
        static_cast<double>(counts[1]) / draw_count
        == Approx(0.3).margin(0.01));
    REQUIRE(
        static_cast<double>(counts[2]) / draw_count
        == Approx(0.5).margin(0.01));
}

TEST_CASE(
    "LogCategoricalDistribution remains stable for extreme log weights",
    "[bayes][stats][log_categorical_distribution]")
{
    using Distribution = LogCategoricalDistribution<3>;
    const Distribution::param_type parameters{Distribution::log_weights_type{
        -std::numeric_limits<double>::infinity(), 1'000.0, 1'002.0}};
    const auto probability_values = parameters.probabilities();
    const Eigen::Map<const Eigen::Vector3d> probabilities{
        probability_values.data()};
    const double denominator = 1.0 + std::exp(-2.0);
    const Eigen::Vector3d expected{
        {0.0, std::exp(-2.0) / denominator, 1.0 / denominator}};

    REQUIRE(probabilities.isApprox(expected));
}

TEST_CASE(
    "LogCategoricalDistribution is invariant to a common log-weight shift",
    "[bayes][stats][log_categorical_distribution]")
{
    using Distribution = LogCategoricalDistribution<3>;
    const Distribution::param_type baseline{
        Distribution::log_weights_type{0.0, 1.0, 2.0}};
    const Distribution::param_type shifted{
        Distribution::log_weights_type{1'000.0, 1'001.0, 1'002.0}};
    const auto baseline_values = baseline.probabilities();
    const auto shifted_values = shifted.probabilities();
    const Eigen::Map<const Eigen::Vector3d> baseline_probabilities{
        baseline_values.data()};
    const Eigen::Map<const Eigen::Vector3d> shifted_probabilities{
        shifted_values.data()};

    REQUIRE(baseline_probabilities.isApprox(shifted_probabilities));
}

TEST_CASE(
    "LogCategoricalDistribution never samples impossible categories",
    "[bayes][stats][log_categorical_distribution]")
{
    using Distribution = LogCategoricalDistribution<3>;
    const Distribution::param_type parameters{Distribution::log_weights_type{
        -std::numeric_limits<double>::infinity(),
        0.0,
        -std::numeric_limits<double>::infinity()}};
    Distribution distribution;
    std::mt19937_64 rng{932};

    for (std::size_t draw_index = 0; draw_index < 100; ++draw_index)
    {
        REQUIRE(distribution(rng, parameters) == 1);
    }
}
