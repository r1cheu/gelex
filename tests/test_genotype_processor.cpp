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

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "gelex/data/genotype/genotype_processor.h"

using namespace gelex;
using namespace gelex::genotype;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace
{
constexpr double k_tolerance = 1e-10;
}  // namespace

TEST_CASE("Standardize<Add> - Basic functionality", "[data]")
{
    SECTION("Happy path - polymorphic variant")
    {
        Eigen::VectorXd variant{{0.0, 1.0, 2.0, 1.0, 0.0}};

        auto stats = Standardize<GeneticMode::A>::process(variant);

        REQUIRE_THAT(stats.mean, WithinRel(0.8, k_tolerance));
        REQUIRE_THAT(stats.stddev, WithinRel(0.8366600265340756, k_tolerance));
        REQUIRE_FALSE(stats.is_monomorphic);

        Eigen::VectorXd expected{
            {-0.9561828874675147,
             0.23904572186687866,
             1.434274331201319,
             0.23904572186687866,
             -0.9561828874675147}};
        REQUIRE(variant.isApprox(expected, k_tolerance));
    }

    SECTION("Happy path - monomorphic variant")
    {
        Eigen::VectorXd variant{{2.0, 2.0, 2.0, 2.0, 2.0}};

        auto stats = Standardize<GeneticMode::A>::process(variant);

        REQUIRE_THAT(stats.mean, WithinRel(2.0, k_tolerance));
        REQUIRE_THAT(stats.stddev, WithinAbs(0.0, k_tolerance));
        REQUIRE(stats.is_monomorphic);
        REQUIRE(variant.isZero(k_tolerance));
    }
}

TEST_CASE("OrthStandardize<Add> - Basic functionality", "[data]")
{
    SECTION("Happy path - polymorphic variant")
    {
        Eigen::VectorXd variant{{0.0, 1.0, 2.0, 1.0, 0.0}};

        auto stats = OrthStandardize<GeneticMode::A>::process(variant);

        REQUIRE_THAT(stats.mean, WithinRel(0.8, k_tolerance));
        double expected_stddev = std::sqrt(0.7);
        REQUIRE_THAT(stats.stddev, WithinRel(expected_stddev, k_tolerance));
        REQUIRE_FALSE(stats.is_monomorphic);

        Eigen::VectorXd expected{
            {-0.9561828874675147,
             0.23904572186687866,
             1.434274331201319,
             0.23904572186687866,
             -0.9561828874675147}};
        REQUIRE(variant.isApprox(expected, k_tolerance));
    }

    SECTION("Happy path - monomorphic variant")
    {
        Eigen::VectorXd variant{{2.0, 2.0, 2.0, 2.0, 2.0}};

        auto stats = OrthStandardize<GeneticMode::A>::process(variant);

        REQUIRE_THAT(stats.mean, WithinRel(2.0, k_tolerance));
        REQUIRE_THAT(stats.stddev, WithinAbs(0.0, k_tolerance));
        REQUIRE(stats.is_monomorphic);
        REQUIRE(variant.isZero(k_tolerance));
    }
}

TEST_CASE("Standardize<Dom> - Basic functionality", "[data]")
{
    SECTION("Happy path - polymorphic variant with heterozygotes")
    {
        Eigen::VectorXd variant{{0.0, 1.0, 2.0, 1.0, 0.0, 2.0}};

        auto stats = Standardize<GeneticMode::D>::process(variant);

        REQUIRE_THAT(stats.mean, WithinRel(0.3333333333333333, k_tolerance));
        REQUIRE_THAT(stats.stddev, WithinRel(0.5163977794943222, k_tolerance));
        REQUIRE_FALSE(stats.is_monomorphic);

        REQUIRE_THAT(variant(2), WithinRel(-0.6454972243679028, k_tolerance));
        REQUIRE_THAT(variant(5), WithinRel(-0.6454972243679028, k_tolerance));
    }

    SECTION("Happy path - variant with no heterozygotes")
    {
        Eigen::VectorXd variant{{0.0, 2.0, 0.0, 2.0}};

        auto stats = Standardize<GeneticMode::D>::process(variant);

        REQUIRE_THAT(stats.mean, WithinRel(0.0, k_tolerance));
        REQUIRE_THAT(stats.stddev, WithinAbs(0.0, k_tolerance));
        REQUIRE(stats.is_monomorphic);
    }
}

TEST_CASE("OrthStandardize<Dom> - Basic functionality", "[data]")
{
    SECTION("Happy path - polymorphic variant")
    {
        Eigen::VectorXd variant{{0.0, 1.0, 2.0, 1.0, 0.0}};

        auto stats = OrthStandardize<GeneticMode::D>::process(variant);

        double expected_mean = 0.24;
        double expected_stddev = std::sqrt(0.288);

        REQUIRE_THAT(stats.mean, WithinRel(expected_mean, k_tolerance));
        REQUIRE_THAT(stats.stddev, WithinRel(expected_stddev, k_tolerance));
        REQUIRE_FALSE(stats.is_monomorphic);

        Eigen::VectorXd expected{
            {-0.24 / expected_stddev,
             0.56 / expected_stddev,
             -0.64 / expected_stddev,
             0.56 / expected_stddev,
             -0.24 / expected_stddev}};
        REQUIRE(variant.isApprox(expected, k_tolerance));
    }

    SECTION("Happy path - monomorphic variant")
    {
        Eigen::VectorXd variant{{2.0, 2.0, 2.0, 2.0, 2.0}};

        auto stats = OrthStandardize<GeneticMode::D>::process(variant);

        REQUIRE_THAT(stats.mean, WithinRel(2.0, k_tolerance));
        REQUIRE_THAT(stats.stddev, WithinAbs(0.0, k_tolerance));
        REQUIRE(stats.is_monomorphic);
        REQUIRE(variant.isZero(k_tolerance));
    }
}

TEST_CASE("StandardizeHWE<Add> uses HWE moments", "[data]")
{
    Eigen::VectorXd variant{{0.0, 1.0, 2.0, 1.0, 0.0}};

    auto stats = StandardizeHWE<GeneticMode::A>::process(variant);

    double expected_mean = 0.8;
    double expected_stddev = std::sqrt(2.0 * 0.4 * 0.6);

    REQUIRE_THAT(stats.mean, WithinRel(expected_mean, k_tolerance));
    REQUIRE_THAT(stats.stddev, WithinRel(expected_stddev, k_tolerance));
    REQUIRE_FALSE(stats.is_monomorphic);

    Eigen::VectorXd expected{
        {-0.8 / expected_stddev,
         0.2 / expected_stddev,
         1.2 / expected_stddev,
         0.2 / expected_stddev,
         -0.8 / expected_stddev}};
    REQUIRE(variant.isApprox(expected, k_tolerance));
}

TEST_CASE("StandardizeHWE<Dom> uses [0,1,0] HWE moments", "[data]")
{
    Eigen::VectorXd variant{{0.0, 1.0, 2.0, 1.0, 0.0}};

    auto stats = StandardizeHWE<GeneticMode::D>::process(variant);

    double expected_mean = 2.0 * 0.4 * 0.6;
    double expected_stddev
        = std::sqrt(2.0 * 0.4 * 0.6 * ((0.4 * 0.4) + (0.6 * 0.6)));

    REQUIRE_THAT(stats.mean, WithinRel(expected_mean, k_tolerance));
    REQUIRE_THAT(stats.stddev, WithinRel(expected_stddev, k_tolerance));
    REQUIRE_FALSE(stats.is_monomorphic);

    Eigen::VectorXd expected{
        {-0.48 / expected_stddev,
         0.52 / expected_stddev,
         -0.48 / expected_stddev,
         0.52 / expected_stddev,
         -0.48 / expected_stddev}};
    REQUIRE(variant.isApprox(expected, k_tolerance));
}

TEST_CASE("OrthStandardizeHWE<Dom> uses [0,2p,4p-2] HWE moments", "[data]")
{
    Eigen::VectorXd variant{{0.0, 1.0, 2.0, 1.0, 0.0}};

    auto stats = OrthStandardizeHWE<GeneticMode::D>::process(variant);

    double expected_mean = 2.0 * 0.4 * 0.4;
    double expected_stddev = 2.0 * 0.4 * 0.6;

    REQUIRE_THAT(stats.mean, WithinRel(expected_mean, k_tolerance));
    REQUIRE_THAT(stats.stddev, WithinRel(expected_stddev, k_tolerance));
    REQUIRE_FALSE(stats.is_monomorphic);

    Eigen::VectorXd expected{
        {-0.32 / expected_stddev,
         0.48 / expected_stddev,
         -0.72 / expected_stddev,
         0.48 / expected_stddev,
         -0.32 / expected_stddev}};
    REQUIRE(variant.isApprox(expected, k_tolerance));
}

TEST_CASE("NOIAStandardize<Add> centers with observed mean", "[data]")
{
    Eigen::VectorXd variant{{0.0, 1.0, 2.0, 1.0, 0.0}};

    auto stats = NOIAStandardize<GeneticMode::A>::process(variant);

    REQUIRE_THAT(stats.mean, WithinRel(0.8, k_tolerance));
    double expected_stddev = std::sqrt(0.7);
    REQUIRE_THAT(stats.stddev, WithinRel(expected_stddev, k_tolerance));
    REQUIRE_FALSE(stats.is_monomorphic);

    Eigen::VectorXd expected{
        {-0.8 / expected_stddev,
         0.2 / expected_stddev,
         1.2 / expected_stddev,
         0.2 / expected_stddev,
         -0.8 / expected_stddev}};
    REQUIRE(variant.isApprox(expected, k_tolerance));
}

TEST_CASE("NOIACenter<Add> centers without scaling", "[data]")
{
    Eigen::VectorXd variant{{0.0, 1.0, 2.0, 1.0, 0.0}};

    auto stats = NOIACenter<GeneticMode::A>::process(variant);

    REQUIRE_THAT(stats.mean, WithinRel(0.8, k_tolerance));
    REQUIRE_FALSE(stats.is_monomorphic);

    Eigen::VectorXd expected{{-0.8, 0.2, 1.2, 0.2, -0.8}};
    REQUIRE(variant.isApprox(expected, k_tolerance));
}

TEST_CASE("NOIACenter<Dom> uses F=4 closed-form with observed freqs", "[data]")
{
    // p_AA=0.2, p_Aa=0.4, p_aa=0.4, D = 0.6 - 0.04 = 0.56
    // dose=2 -> -2*p_aa*p_Aa/D = -0.32/0.56
    // dose=1 ->  4*p_AA*p_aa/D =  0.32/0.56
    // dose=0 -> -2*p_AA*p_Aa/D = -0.16/0.56
    Eigen::VectorXd variant{{0.0, 1.0, 2.0, 1.0, 0.0}};

    auto stats = NOIACenter<GeneticMode::D>::process(variant);

    double c_AA = -0.32 / 0.56;
    double c_Aa = 0.32 / 0.56;
    double c_aa = -0.16 / 0.56;

    REQUIRE_THAT(stats.mean, WithinAbs(0.0, k_tolerance));
    REQUIRE_FALSE(stats.is_monomorphic);

    Eigen::VectorXd expected{{c_aa, c_Aa, c_AA, c_Aa, c_aa}};
    REQUIRE(variant.isApprox(expected, k_tolerance));
}

TEST_CASE("NOIAStandardize<Dom> scales to unit variance", "[data]")
{
    Eigen::VectorXd variant{{0.0, 1.0, 2.0, 1.0, 0.0}};

    auto stats = NOIAStandardize<GeneticMode::D>::process(variant);

    double c_AA = -0.32 / 0.56;
    double c_Aa = 0.32 / 0.56;
    double c_aa = -0.16 / 0.56;
    double sum_sq = (2.0 * c_aa * c_aa) + (2.0 * c_Aa * c_Aa) + (c_AA * c_AA);
    double expected_stddev = std::sqrt(sum_sq / 4.0);

    REQUIRE_THAT(stats.stddev, WithinRel(expected_stddev, k_tolerance));
    REQUIRE_FALSE(stats.is_monomorphic);

    Eigen::VectorXd expected{
        {c_aa / expected_stddev,
         c_Aa / expected_stddev,
         c_AA / expected_stddev,
         c_Aa / expected_stddev,
         c_aa / expected_stddev}};
    REQUIRE(variant.isApprox(expected, k_tolerance));
}

TEST_CASE("NOIA additive and dominance columns are sample-orthogonal", "[data]")
{
    // Sample deviating from HWE: p_AA=0.3, p_Aa=0.1, p_aa=0.6 (MAF=0.35)
    Eigen::VectorXd add_variant{
        {2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    Eigen::VectorXd dom_variant = add_variant;

    NOIACenter<GeneticMode::A>::process(add_variant);
    NOIACenter<GeneticMode::D>::process(dom_variant);

    double inner = add_variant.dot(dom_variant);
    REQUIRE_THAT(inner, WithinAbs(0.0, k_tolerance));
}

TEST_CASE("NOIA flags monomorphic samples", "[data]")
{
    SECTION("All homozygous-alt")
    {
        Eigen::VectorXd variant{{2.0, 2.0, 2.0, 2.0, 2.0}};
        auto stats = NOIAStandardize<GeneticMode::D>::process(variant);
        REQUIRE(stats.is_monomorphic);
    }

    SECTION("No heterozygotes, both homozygous classes present")
    {
        Eigen::VectorXd variant{{0.0, 2.0, 0.0, 2.0}};
        auto stats = NOIAStandardize<GeneticMode::D>::process(variant);
        // p_Aa=0 makes every NOIA dominance coefficient zero.
        REQUIRE(stats.is_monomorphic);
    }

    SECTION("Missing AA class")
    {
        Eigen::VectorXd variant{{0.0, 1.0, 0.0, 1.0, 0.0}};
        auto stats = NOIAStandardize<GeneticMode::D>::process(variant);
        REQUIRE(stats.is_monomorphic);
    }
}

TEST_CASE("NOIA handles NaN by imputing to analytic mean", "[data]")
{
    double nan_v = std::nan("");
    Eigen::VectorXd variant{{0.0, 1.0, nan_v, 2.0, 1.0}};

    auto stats = NOIAStandardize<GeneticMode::A>::process(variant);

    REQUIRE_FALSE(stats.is_monomorphic);
    // NaN was imputed to mean, then centered -> 0.
    REQUIRE_THAT(variant(2), WithinAbs(0.0, k_tolerance));
}
