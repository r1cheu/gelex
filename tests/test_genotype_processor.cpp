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

#include <cmath>
#include <initializer_list>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "gelex/data/genotype_processor.h"

using namespace gelex;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace
{

constexpr double k_tolerance = 1e-10;

auto make_vector(std::initializer_list<double> values) -> Eigen::VectorXd
{
    Eigen::VectorXd vector(static_cast<Eigen::Index>(values.size()));
    Eigen::Index index = 0;
    for (double value : values)
    {
        vector(index) = value;
        ++index;
    }
    return vector;
}

auto require_vector_within_rel(
    Eigen::Ref<const Eigen::VectorXd> actual,
    Eigen::Ref<const Eigen::VectorXd> expected,
    double tolerance = k_tolerance) -> void
{
    REQUIRE(actual.size() == expected.size());
    for (Eigen::Index i = 0; i < actual.size(); ++i)
    {
        REQUIRE_THAT(actual(i), WithinRel(expected(i), tolerance));
    }
}

auto require_vector_within_abs(
    Eigen::Ref<const Eigen::VectorXd> actual,
    double expected,
    double tolerance = k_tolerance) -> void
{
    for (Eigen::Index i = 0; i < actual.size(); ++i)
    {
        REQUIRE_THAT(actual(i), WithinAbs(expected, tolerance));
    }
}

}  // namespace

TEST_CASE(
    "AdditiveProcessor<StandardizeMethod> - Basic functionality",
    "[data]")
{
    SECTION("Happy path - polymorphic variant")
    {
        Eigen::VectorXd variant = make_vector({0.0, 1.0, 2.0, 1.0, 0.0});

        auto stats
            = AdditiveProcessor<StandardizeMethod>::process_variant(variant);

        REQUIRE_THAT(stats.mean, WithinRel(0.8, k_tolerance));
        REQUIRE_THAT(stats.stddev, WithinRel(0.8366600265340756, k_tolerance));
        REQUIRE_FALSE(stats.is_monomorphic);

        Eigen::VectorXd expected = make_vector(
            {-0.9561828874675147,
             0.23904572186687866,
             1.434274331201319,
             0.23904572186687866,
             -0.9561828874675147});
        require_vector_within_rel(variant, expected);
    }

    SECTION("Happy path - monomorphic variant")
    {
        Eigen::VectorXd variant = make_vector({2.0, 2.0, 2.0, 2.0, 2.0});

        auto stats
            = AdditiveProcessor<StandardizeMethod>::process_variant(variant);

        REQUIRE_THAT(stats.mean, WithinRel(2.0, k_tolerance));
        REQUIRE_THAT(stats.stddev, WithinAbs(0.0, k_tolerance));
        REQUIRE(stats.is_monomorphic);
        require_vector_within_abs(variant, 0.0);
    }
}

TEST_CASE(
    "AdditiveProcessor<OrthStandardizeMethod> - Basic functionality",
    "[data]")
{
    SECTION("Happy path - polymorphic variant")
    {
        Eigen::VectorXd variant = make_vector({0.0, 1.0, 2.0, 1.0, 0.0});

        auto stats = AdditiveProcessor<OrthStandardizeMethod>::process_variant(
            variant);

        REQUIRE_THAT(stats.mean, WithinRel(0.8, k_tolerance));
        double expected_stddev = std::sqrt(0.7);
        REQUIRE_THAT(stats.stddev, WithinRel(expected_stddev, k_tolerance));
        REQUIRE_FALSE(stats.is_monomorphic);

        Eigen::VectorXd expected = make_vector(
            {-0.9561828874675147,
             0.23904572186687866,
             1.434274331201319,
             0.23904572186687866,
             -0.9561828874675147});
        require_vector_within_rel(variant, expected);
    }

    SECTION("Happy path - monomorphic variant")
    {
        Eigen::VectorXd variant = make_vector({2.0, 2.0, 2.0, 2.0, 2.0});

        auto stats = AdditiveProcessor<OrthStandardizeMethod>::process_variant(
            variant);

        REQUIRE_THAT(stats.mean, WithinRel(2.0, k_tolerance));
        REQUIRE_THAT(stats.stddev, WithinAbs(0.0, k_tolerance));
        REQUIRE(stats.is_monomorphic);
        require_vector_within_abs(variant, 0.0);
    }
}

TEST_CASE(
    "DominantProcessor<StandardizeMethod> - Basic functionality",
    "[data]")
{
    SECTION("Happy path - polymorphic variant with heterozygotes")
    {
        Eigen::VectorXd variant = make_vector({0.0, 1.0, 2.0, 1.0, 0.0, 2.0});

        auto stats
            = DominantProcessor<StandardizeMethod>::process_variant(variant);

        REQUIRE_THAT(stats.mean, WithinRel(0.3333333333333333, k_tolerance));
        REQUIRE_THAT(stats.stddev, WithinRel(0.5163977794943222, k_tolerance));
        REQUIRE_FALSE(stats.is_monomorphic);

        REQUIRE_THAT(variant(2), WithinRel(-0.6454972243679028, k_tolerance));
        REQUIRE_THAT(variant(5), WithinRel(-0.6454972243679028, k_tolerance));
    }

    SECTION("Happy path - variant with no heterozygotes")
    {
        Eigen::VectorXd variant = make_vector({0.0, 2.0, 0.0, 2.0});

        auto stats
            = DominantProcessor<StandardizeMethod>::process_variant(variant);

        REQUIRE_THAT(stats.mean, WithinRel(0.0, k_tolerance));
        REQUIRE_THAT(stats.stddev, WithinAbs(0.0, k_tolerance));
        REQUIRE(stats.is_monomorphic);
    }
}

TEST_CASE(
    "DominantProcessor<OrthStandardizeMethod> - Basic functionality",
    "[data]")
{
    SECTION("Happy path - polymorphic variant")
    {
        Eigen::VectorXd variant = make_vector({0.0, 1.0, 2.0, 1.0, 0.0});

        auto stats = DominantProcessor<OrthStandardizeMethod>::process_variant(
            variant);

        double expected_mean = 0.24;
        double expected_stddev = std::sqrt(0.288);

        REQUIRE_THAT(stats.mean, WithinRel(expected_mean, k_tolerance));
        REQUIRE_THAT(stats.stddev, WithinRel(expected_stddev, k_tolerance));
        REQUIRE_FALSE(stats.is_monomorphic);

        Eigen::VectorXd expected = make_vector(
            {-0.24 / expected_stddev,
             0.56 / expected_stddev,
             -0.64 / expected_stddev,
             0.56 / expected_stddev,
             -0.24 / expected_stddev});
        require_vector_within_rel(variant, expected);
    }

    SECTION("Happy path - monomorphic variant")
    {
        Eigen::VectorXd variant = make_vector({2.0, 2.0, 2.0, 2.0, 2.0});

        auto stats = DominantProcessor<OrthStandardizeMethod>::process_variant(
            variant);

        REQUIRE_THAT(stats.mean, WithinRel(2.0, k_tolerance));
        REQUIRE_THAT(stats.stddev, WithinAbs(0.0, k_tolerance));
        REQUIRE(stats.is_monomorphic);
        require_vector_within_abs(variant, 0.0);
    }
}

TEST_CASE("AdditiveProcessor<StandardizeHWEMethod> uses HWE moments", "[data]")
{
    Eigen::VectorXd variant = make_vector({0.0, 1.0, 2.0, 1.0, 0.0});

    auto stats
        = AdditiveProcessor<StandardizeHWEMethod>::process_variant(variant);

    double expected_mean = 0.8;
    double expected_stddev = std::sqrt(2.0 * 0.4 * 0.6);

    REQUIRE_THAT(stats.mean, WithinRel(expected_mean, k_tolerance));
    REQUIRE_THAT(stats.stddev, WithinRel(expected_stddev, k_tolerance));
    REQUIRE_FALSE(stats.is_monomorphic);

    Eigen::VectorXd expected = make_vector(
        {-0.8 / expected_stddev,
         0.2 / expected_stddev,
         1.2 / expected_stddev,
         0.2 / expected_stddev,
         -0.8 / expected_stddev});
    require_vector_within_rel(variant, expected);
}

TEST_CASE(
    "DominantProcessor<StandardizeHWEMethod> uses [0,1,0] HWE moments",
    "[data]")
{
    Eigen::VectorXd variant = make_vector({0.0, 1.0, 2.0, 1.0, 0.0});

    auto stats
        = DominantProcessor<StandardizeHWEMethod>::process_variant(variant);

    double expected_mean = 2.0 * 0.4 * 0.6;
    double expected_stddev
        = std::sqrt(2.0 * 0.4 * 0.6 * ((0.4 * 0.4) + (0.6 * 0.6)));

    REQUIRE_THAT(stats.mean, WithinRel(expected_mean, k_tolerance));
    REQUIRE_THAT(stats.stddev, WithinRel(expected_stddev, k_tolerance));
    REQUIRE_FALSE(stats.is_monomorphic);

    Eigen::VectorXd expected = make_vector(
        {-0.48 / expected_stddev,
         0.52 / expected_stddev,
         -0.48 / expected_stddev,
         0.52 / expected_stddev,
         -0.48 / expected_stddev});
    require_vector_within_rel(variant, expected);
}

TEST_CASE(
    "DominantProcessor<OrthStandardizeHWEMethod> uses [0,2p,4p-2] HWE moments",
    "[data]")
{
    Eigen::VectorXd variant = make_vector({0.0, 1.0, 2.0, 1.0, 0.0});

    auto stats
        = DominantProcessor<OrthStandardizeHWEMethod>::process_variant(variant);

    double expected_mean = 2.0 * 0.4 * 0.4;
    double expected_stddev = 2.0 * 0.4 * 0.6;

    REQUIRE_THAT(stats.mean, WithinRel(expected_mean, k_tolerance));
    REQUIRE_THAT(stats.stddev, WithinRel(expected_stddev, k_tolerance));
    REQUIRE_FALSE(stats.is_monomorphic);

    Eigen::VectorXd expected = make_vector(
        {-0.32 / expected_stddev,
         0.48 / expected_stddev,
         -0.72 / expected_stddev,
         0.48 / expected_stddev,
         -0.32 / expected_stddev});
    require_vector_within_rel(variant, expected);
}
