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

#include "gelex/data/genotype_processor.h"

using namespace gelex;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("AdditiveProcessor<StandardizeMethod> - Basic functionality", "[data]")
{
    SECTION("Happy path - polymorphic variant")
    {
        Eigen::VectorXd variant(5);
        variant << 0.0, 1.0, 2.0, 1.0, 0.0;
        Eigen::VectorXd original = variant;

        auto stats = AdditiveProcessor<StandardizeMethod>::process_variant(variant);

        // Check statistics
        REQUIRE_THAT(stats.mean, WithinRel(0.8, 1e-10));
        REQUIRE_THAT(stats.stddev, WithinRel(0.8366600265340756, 1e-10));
        REQUIRE_FALSE(stats.is_monomorphic);

        // Check standardization
        Eigen::VectorXd expected(5);
        expected << -0.9561828874675147, 0.23904572186687866, 1.434274331201319,
            0.23904572186687866, -0.9561828874675147;

        for (int i = 0; i < variant.size(); ++i)
        {
            REQUIRE_THAT(variant(i), WithinRel(expected(i), 1e-10));
        }
    }

    SECTION("Happy path - monomorphic variant")
    {
        Eigen::VectorXd variant(5);
        variant << 2.0, 2.0, 2.0, 2.0, 2.0;

        auto stats = AdditiveProcessor<StandardizeMethod>::process_variant(variant);

        REQUIRE_THAT(stats.mean, WithinRel(2.0, 1e-10));
        REQUIRE_THAT(stats.stddev, WithinAbs(0.0, 1e-10));
        REQUIRE(stats.is_monomorphic);

        // CenterMethod centers to [0,0,0,0,0], monomorphic so not divided
        for (int i = 0; i < variant.size(); ++i)
        {
            REQUIRE_THAT(variant(i), WithinAbs(0.0, 1e-10));
        }
    }
}

TEST_CASE("AdditiveProcessor<OrthStandardizeMethod> - Basic functionality", "[data]")
{
    SECTION("Happy path - polymorphic variant")
    {
        Eigen::VectorXd variant(5);
        variant << 0.0, 1.0, 2.0, 1.0, 0.0;
        Eigen::VectorXd original = variant;

        auto stats = AdditiveProcessor<OrthStandardizeMethod>::process_variant(variant);

        // OrthStandardizeMethod = CenterMethod + divide by sample stddev
        // Same as StandardizeMethod for additive
        REQUIRE_THAT(stats.mean, WithinRel(0.8, 1e-10));
        double expected_stddev = std::sqrt(0.7);
        REQUIRE_THAT(stats.stddev, WithinRel(expected_stddev, 1e-10));
        REQUIRE_FALSE(stats.is_monomorphic);

        // values = (x - mean) / stddev = same as StandardizeMethod
        Eigen::VectorXd expected(5);
        expected << -0.9561828874675147, 0.23904572186687866,
            1.434274331201319, 0.23904572186687866, -0.9561828874675147;

        for (int i = 0; i < variant.size(); ++i)
        {
            REQUIRE_THAT(variant(i), WithinRel(expected(i), 1e-10));
        }
    }

    SECTION("Happy path - monomorphic variant")
    {
        Eigen::VectorXd variant(5);
        variant << 2.0, 2.0, 2.0, 2.0, 2.0;

        auto stats = AdditiveProcessor<OrthStandardizeMethod>::process_variant(variant);

        REQUIRE_THAT(stats.mean, WithinRel(2.0, 1e-10));
        REQUIRE_THAT(stats.stddev, WithinAbs(0.0, 1e-10));
        REQUIRE(stats.is_monomorphic);

        // CenterMethod centers to [0,0,...], monomorphic so not divided
        for (int i = 0; i < variant.size(); ++i)
        {
            REQUIRE_THAT(variant(i), WithinAbs(0.0, 1e-10));
        }
    }
}

TEST_CASE("DominantProcessor<StandardizeMethod> - Basic functionality", "[data]")
{
    SECTION("Happy path - polymorphic variant with heterozygotes")
    {
        Eigen::VectorXd variant(6);
        variant << 0.0, 1.0, 2.0, 1.0, 0.0, 2.0;

        auto stats = DominantProcessor<StandardizeMethod>::process_variant(variant);

        // After converting 2.0 to 0.0, we have: 0.0, 1.0, 0.0, 1.0, 0.0, 0.0
        // Mean = (0+1+0+1+0+0)/6 = 2/6 = 0.333333
        // Variance = [(0-0.333)^2*4 + (1-0.333)^2*2]/5 = [0.111111*4 +
        // 0.444444*2]/5 = [0.444444 + 0.888888]/5 = 1.333332/5 = 0.2666664
        // Stddev = sqrt(0.2666664) = 0.5163978

        REQUIRE_THAT(stats.mean, WithinRel(0.3333333333333333, 1e-10));
        REQUIRE_THAT(stats.stddev, WithinRel(0.5163977794943222, 1e-10));
        REQUIRE_FALSE(stats.is_monomorphic);

        // Check that 2.0 values were converted to 0.0 and then standardized
        // Original value 2.0 -> 0.0, then (0.0 - 0.333333)/0.5163978 =
        // -0.645497
        REQUIRE_THAT(variant(2), WithinRel(-0.6454972243679028, 1e-10));
        REQUIRE_THAT(variant(5), WithinRel(-0.6454972243679028, 1e-10));
    }

    SECTION("Happy path - variant with no heterozygotes")
    {
        Eigen::VectorXd variant(4);
        variant << 0.0, 2.0, 0.0, 2.0;

        auto stats = DominantProcessor<StandardizeMethod>::process_variant(variant);

        // After converting 2.0 to 0.0, we have all zeros
        REQUIRE_THAT(stats.mean, WithinRel(0.0, 1e-10));
        REQUIRE_THAT(stats.stddev, WithinAbs(0.0, 1e-10));
        REQUIRE(stats.is_monomorphic);
    }
}

TEST_CASE("DominantProcessor<OrthStandardizeMethod> - Basic functionality", "[data]")
{
    SECTION("Happy path - polymorphic variant")
    {
        Eigen::VectorXd variant(5);
        variant << 0.0, 1.0, 2.0, 1.0, 0.0;

        auto stats = DominantProcessor<OrthStandardizeMethod>::process_variant(variant);

        // OrthCenterMethod recodes: maf=0.4
        // 0→0.0, 1→2*0.4=0.8, 2→4*0.4-2=-0.4
        // recoded: [0, 0.8, -0.4, 0.8, 0]
        // CenterMethod: mean = (0+0.8-0.4+0.8+0)/5 = 1.2/5 = 0.24
        // centered: [-0.24, 0.56, -0.64, 0.56, -0.24]
        // sum_sq = 0.0576+0.3136+0.4096+0.3136+0.0576 = 1.152
        // sample_var = 1.152/4 = 0.288, stddev = sqrt(0.288)

        double expected_mean = 0.24;
        double expected_stddev = std::sqrt(0.288);

        REQUIRE_THAT(stats.mean, WithinRel(expected_mean, 1e-10));
        REQUIRE_THAT(stats.stddev, WithinRel(expected_stddev, 1e-10));
        REQUIRE_FALSE(stats.is_monomorphic);

        // divided by stddev
        Eigen::VectorXd expected(5);
        expected << -0.24 / expected_stddev, 0.56 / expected_stddev,
            -0.64 / expected_stddev, 0.56 / expected_stddev,
            -0.24 / expected_stddev;

        for (int i = 0; i < variant.size(); ++i)
        {
            REQUIRE_THAT(variant(i), WithinRel(expected(i), 1e-10));
        }
    }

    SECTION("Happy path - monomorphic variant")
    {
        Eigen::VectorXd variant(5);
        variant << 2.0, 2.0, 2.0, 2.0, 2.0;

        auto stats = DominantProcessor<OrthStandardizeMethod>::process_variant(variant);

        // maf = 2.0/2 = 1.0
        // two_alt_encode = 4*1.0 - 2 = 2.0
        // recoded: [2, 2, 2, 2, 2]
        // CenterMethod centers: mean=2.0, centered to [0,0,0,0,0]
        // stddev=0, is_monomorphic=true

        REQUIRE_THAT(stats.mean, WithinRel(2.0, 1e-10));
        REQUIRE_THAT(stats.stddev, WithinAbs(0.0, 1e-10));
        REQUIRE(stats.is_monomorphic);

        for (int i = 0; i < variant.size(); ++i)
        {
            REQUIRE_THAT(variant(i), WithinAbs(0.0, 1e-10));
        }
    }
}
