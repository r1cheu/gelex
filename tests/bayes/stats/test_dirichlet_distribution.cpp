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
#include <cstddef>
#include <random>
#include <type_traits>

#include "gelex/bayes/stats/dirichlet_distribution.h"

using Catch::Approx;
using gelex::DirichletDistribution;

template <std::size_t K>
concept ValidDirichletDistribution
    = requires { typename DirichletDistribution<K>; };

static_assert(ValidDirichletDistribution<2>);
static_assert(!ValidDirichletDistribution<1>);

TEST_CASE(
    "DirichletDistribution follows the parameter API",
    "[bayes][stats][dirichlet_distribution]")
{
    using Distribution = DirichletDistribution<3>;
    using Parameters = Distribution::param_type;
    static_assert(std::is_same_v<Parameters::distribution_type, Distribution>);

    const Parameters default_parameters;
    const Distribution::result_type expected_default{1.0, 1.0, 1.0};
    REQUIRE(default_parameters.concentrations() == expected_default);

    const Parameters parameters{Distribution::result_type{1.0, 2.0, 3.0}};
    Distribution distribution;
    distribution.param(parameters);
    REQUIRE(distribution.param() == parameters);
    REQUIRE(distribution.concentrations() == parameters.concentrations());
    REQUIRE(Distribution{parameters}.param() == parameters);
    REQUIRE(Distribution{parameters.concentrations()}.param() == parameters);

    distribution.reset();
    REQUIRE(distribution.param() == parameters);

    std::mt19937_64 stored_rng{481};
    std::mt19937_64 supplied_rng{481};
    Distribution supplied_distribution;
    for (std::size_t draw_index = 0; draw_index < 100; ++draw_index)
    {
        const auto stored_sample = distribution(stored_rng);
        const auto supplied_sample
            = supplied_distribution(supplied_rng, parameters);
        const Eigen::Map<const Eigen::Vector3d> stored_values{
            stored_sample.data()};
        const Eigen::Map<const Eigen::Vector3d> supplied_values{
            supplied_sample.data()};
        REQUIRE(stored_values.isApprox(supplied_values));
    }
}

TEST_CASE(
    "DirichletDistribution samples the requested concentrations",
    "[bayes][stats][dirichlet_distribution]")
{
    constexpr std::size_t draw_count = 50'000;
    using Distribution = DirichletDistribution<3>;
    const Distribution::param_type parameters{
        Distribution::result_type{1.0, 2.0, 3.0}};
    Distribution distribution;
    std::mt19937_64 rng{123};
    Eigen::Vector3d sum = Eigen::Vector3d::Zero();
    Eigen::Matrix3d sum_outer_products = Eigen::Matrix3d::Zero();

    for (std::size_t draw_index = 0; draw_index < draw_count; ++draw_index)
    {
        const auto sample = distribution(rng, parameters);
        const Eigen::Map<const Eigen::Vector3d> values{sample.data()};
        sum += values;
        sum_outer_products.noalias() += values * values.transpose();
    }

    const auto count = static_cast<double>(draw_count);
    const Eigen::Vector3d mean = sum / count;
    const Eigen::Matrix3d covariance
        = (sum_outer_products - (count * mean * mean.transpose()))
          / (count - 1.0);
    const Eigen::Vector3d expected_mean{{1.0 / 6.0, 2.0 / 6.0, 3.0 / 6.0}};
    const Eigen::Matrix3d expected_covariance{
        {5.0 / 252.0, -2.0 / 252.0, -3.0 / 252.0},
        {-2.0 / 252.0, 8.0 / 252.0, -6.0 / 252.0},
        {-3.0 / 252.0, -6.0 / 252.0, 9.0 / 252.0}};

    REQUIRE(mean.isApprox(expected_mean, 0.01));
    REQUIRE(covariance.isApprox(expected_covariance, 0.03));
    REQUIRE(mean.sum() == Approx(1.0));
}
