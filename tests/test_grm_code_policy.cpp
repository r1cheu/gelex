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

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/data/genotype/genotype_processor.h"

using namespace gelex;

// ============================================================================
// AdditiveProcessor<CenterMethod>
// ============================================================================

TEST_CASE(
    "AdditiveProcessor<CenterMethod> single column via process_matrix",
    "[grm]")
{
    Eigen::MatrixXd geno(5, 1);
    geno << 0, 1, 2, 1, 0;
    // mean = (0+1+2+1+0)/5 = 0.8

    Eigen::MatrixXd expected(5, 1);
    expected << -0.8, 0.2, 1.2, 0.2, -0.8;

    process_matrix<GeneticEffectType::Add>(GenotypeProcessMethod::Center, geno);

    REQUIRE(geno.isApprox(expected));
}

TEST_CASE(
    "AdditiveProcessor<CenterMethod> multiple columns via process_matrix",
    "[grm]")
{
    Eigen::MatrixXd geno(4, 3);
    geno << 0, 2, 1, 1, 1, 1, 2, 0, 1, 1, 1, 1;
    // col0: mean = 1.0
    // col1: mean = 1.0
    // col2: mean = 1.0

    Eigen::MatrixXd expected(4, 3);
    expected << -1, 1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0;

    process_matrix<GeneticEffectType::Add>(GenotypeProcessMethod::Center, geno);

    REQUIRE(geno.isApprox(expected));
}

TEST_CASE(
    "AdditiveProcessor<CenterMethod> already centered via process_matrix",
    "[grm]")
{
    Eigen::MatrixXd geno(4, 1);
    geno << -1, 1, -1, 1;  // mean = 0

    Eigen::MatrixXd expected = geno;

    process_matrix<GeneticEffectType::Add>(GenotypeProcessMethod::Center, geno);

    REQUIRE(geno.isApprox(expected));
}

// ============================================================================
// Center additive (= AdditiveProcessor<CenterMethod> via process_matrix)
// ============================================================================

TEST_CASE("Center additive mode", "[grm][Center]")
{
    SECTION("basic centering")
    {
        Eigen::MatrixXd geno(5, 1);
        geno << 0, 1, 2, 1, 0;
        // mean = 0.8

        Eigen::MatrixXd expected(5, 1);
        expected << -0.8, 0.2, 1.2, 0.2, -0.8;

        process_matrix<GeneticEffectType::Add>(
            GenotypeProcessMethod::Center, geno);

        REQUIRE(geno.isApprox(expected));
    }

    SECTION("multiple columns")
    {
        Eigen::MatrixXd geno(4, 2);
        geno << 0, 2, 1, 1, 2, 0, 1, 1;
        // col0: mean=1, col1: mean=1

        Eigen::MatrixXd expected(4, 2);
        expected << -1, 1, 0, 0, 1, -1, 0, 0;

        process_matrix<GeneticEffectType::Add>(
            GenotypeProcessMethod::Center, geno);

        REQUIRE(geno.isApprox(expected));
    }
}

// ============================================================================
// Center dominance (= DominantProcessor<CenterMethod> via process_matrix)
// ============================================================================

TEST_CASE("Center dominance mode", "[grm][Center]")
{
    SECTION("basic transformation")
    {
        Eigen::MatrixXd geno(5, 1);
        geno << 0, 1, 2, 1, 0;

        // 2→0: [0, 1, 0, 1, 0], mean=0.4, centered: [-0.4, 0.6, -0.4, 0.6,
        // -0.4]
        Eigen::MatrixXd expected(5, 1);
        expected << -0.4, 0.6, -0.4, 0.6, -0.4;

        process_matrix<GeneticEffectType::Dom>(
            GenotypeProcessMethod::Center, geno);

        REQUIRE(geno.isApprox(expected, 1e-10));
    }

    SECTION("all heterozygous")
    {
        // all 1s: 2→0 gives [1,1,1,1], mean=1.0, centered: [0,0,0,0]
        Eigen::MatrixXd geno(4, 1);
        geno << 1, 1, 1, 1;

        Eigen::MatrixXd expected(4, 1);
        expected << 0, 0, 0, 0;

        process_matrix<GeneticEffectType::Dom>(
            GenotypeProcessMethod::Center, geno);

        REQUIRE(geno.isApprox(expected, 1e-10));
    }
}

// ============================================================================
// OrthCenter additive (= AdditiveProcessor<OrthCenterMethod> via
// process_matrix)
// ============================================================================

TEST_CASE("OrthCenter additive mode", "[grm][OrthCenter]")
{
    SECTION("basic centering")
    {
        Eigen::MatrixXd geno(5, 1);
        geno << 0, 1, 2, 1, 0;

        Eigen::MatrixXd expected(5, 1);
        expected << -0.8, 0.2, 1.2, 0.2, -0.8;

        process_matrix<GeneticEffectType::Add>(
            GenotypeProcessMethod::OrthCenter, geno);

        REQUIRE(geno.isApprox(expected));
    }
}

// ============================================================================
// OrthCenter dominance (= DominantProcessor<OrthCenterMethod> via
// process_matrix)
// ============================================================================

TEST_CASE("OrthCenter dominance mode", "[grm][OrthCenter]")
{
    SECTION("basic transformation")
    {
        // genotype: 0, 1, 2, 1, 0
        // maf = 0.8/2 = 0.4
        // orthogonal Dom recodes: 0→0.0, 1→2*0.4=0.8, 2→4*0.4-2=-0.4
        // recoded: [0, 0.8, -0.4, 0.8, 0]
        // sample mean = (0+0.8-0.4+0.8+0)/5 = 0.24
        // centered: [-0.24, 0.56, -0.64, 0.56, -0.24]
        Eigen::MatrixXd geno(5, 1);
        geno << 0, 1, 2, 1, 0;

        Eigen::MatrixXd expected(5, 1);
        expected << -0.24, 0.56, -0.64, 0.56, -0.24;

        process_matrix<GeneticEffectType::Dom>(
            GenotypeProcessMethod::OrthCenter, geno);

        REQUIRE(geno.isApprox(expected, 1e-10));
    }

    SECTION("all heterozygous")
    {
        // all 1s: maf = 0.5
        // OrthCenter recodes: 1→2*0.5=1.0
        // recoded: [1, 1, 1, 1]
        // CenterMethod: mean=1.0, centered: [0, 0, 0, 0]
        Eigen::MatrixXd geno(4, 1);
        geno << 1, 1, 1, 1;

        Eigen::MatrixXd expected(4, 1);
        expected << 0, 0, 0, 0;

        process_matrix<GeneticEffectType::Dom>(
            GenotypeProcessMethod::OrthCenter, geno);

        REQUIRE(geno.isApprox(expected, 1e-10));
    }

    SECTION("all homozygous AA")
    {
        // all 2s: maf = 1.0
        // OrthCenter recodes: 2→4*1.0-2=2.0
        // recoded: [2, 2, 2, 2]
        // CenterMethod: mean=2.0, centered: [0, 0, 0, 0]
        Eigen::MatrixXd geno(4, 1);
        geno << 2, 2, 2, 2;

        Eigen::MatrixXd expected(4, 1);
        expected << 0, 0, 0, 0;

        process_matrix<GeneticEffectType::Dom>(
            GenotypeProcessMethod::OrthCenter, geno);

        REQUIRE(geno.isApprox(expected));
    }
}

// ============================================================================
// OrthCenterHWE dominance
// ============================================================================

TEST_CASE("OrthCenterHWE dominance mode", "[grm][OrthCenterHWE]")
{
    SECTION("basic transformation")
    {
        // genotype: 0, 1, 2, 1, 0
        // maf = 0.8/2 = 0.4
        // orthogonal Dom recodes: 0→0.0, 1→2*0.4=0.8, 2→4*0.4-2=-0.4
        // recoded: [0, 0.8, -0.4, 0.8, 0]
        // OrthHWE mean = 2*maf² = 2*0.16 = 0.32
        // centered: [-0.32, 0.48, -0.72, 0.48, -0.32]
        Eigen::MatrixXd geno(5, 1);
        geno << 0, 1, 2, 1, 0;

        Eigen::MatrixXd expected(5, 1);
        expected << -0.32, 0.48, -0.72, 0.48, -0.32;

        process_matrix<GeneticEffectType::Dom>(
            GenotypeProcessMethod::OrthCenterHWE, geno);

        REQUIRE(geno.isApprox(expected, 1e-10));
    }
}

// ============================================================================
// OrthStandardize additive (= AdditiveProcessor<OrthStandardizeMethod> via
// process_matrix)
// ============================================================================

TEST_CASE("OrthStandardize additive mode", "[grm][OrthStandardize]")
{
    SECTION("basic standardization")
    {
        // genotype: 0, 1, 2, 1, 0
        // CenterMethod: mean=0.8, centered=[-0.8, 0.2, 1.2, 0.2, -0.8]
        // sample stddev = sqrt(sum_sq/4) = sqrt(2.8/4) = sqrt(0.7)
        // divided by stddev
        Eigen::MatrixXd geno(5, 1);
        geno << 0, 1, 2, 1, 0;

        double mean = 0.8;
        double stddev = std::sqrt(0.7);

        Eigen::MatrixXd expected(5, 1);
        expected << (0 - mean) / stddev, (1 - mean) / stddev,
            (2 - mean) / stddev, (1 - mean) / stddev, (0 - mean) / stddev;

        process_matrix<GeneticEffectType::Add>(
            GenotypeProcessMethod::OrthStandardize, geno);

        REQUIRE(geno.isApprox(expected));
    }

    SECTION("monomorphic SNP sets to zero")
    {
        Eigen::MatrixXd geno(4, 1);
        geno << 0, 0, 0, 0;

        process_matrix<GeneticEffectType::Add>(
            GenotypeProcessMethod::OrthStandardize, geno);

        REQUIRE(geno.isApproxToConstant(0.0));
    }

    SECTION("all homozygous AA sets to zero")
    {
        Eigen::MatrixXd geno(4, 1);
        geno << 2, 2, 2, 2;

        process_matrix<GeneticEffectType::Add>(
            GenotypeProcessMethod::OrthStandardize, geno);

        REQUIRE(geno.isApproxToConstant(0.0));
    }
}

// ============================================================================
// OrthStandardize dominance (= DominantProcessor<OrthStandardizeMethod> via
// process_matrix)
// ============================================================================

TEST_CASE("OrthStandardize dominance mode", "[grm][OrthStandardize]")
{
    SECTION("basic transformation")
    {
        // genotype: 0, 1, 2, 1, 0
        // maf=0.4, orthogonal Dom recodes: 0→0, 1→0.8, 2→-0.4
        // recoded: [0, 0.8, -0.4, 0.8, 0]
        // sample mean=0.24, centered: [-0.24, 0.56, -0.64, 0.56, -0.24]
        // sample stddev = sqrt(1.152/4) = sqrt(0.288)
        // divided by stddev
        Eigen::MatrixXd geno(5, 1);
        geno << 0, 1, 2, 1, 0;

        double stddev = std::sqrt(0.288);

        Eigen::MatrixXd expected(5, 1);
        expected << -0.24 / stddev, 0.56 / stddev, -0.64 / stddev,
            0.56 / stddev, -0.24 / stddev;

        process_matrix<GeneticEffectType::Dom>(
            GenotypeProcessMethod::OrthStandardize, geno);

        REQUIRE(geno.isApprox(expected));
    }

    SECTION("monomorphic SNP sets to zero")
    {
        Eigen::MatrixXd geno(4, 1);
        geno << 0, 0, 0, 0;

        process_matrix<GeneticEffectType::Dom>(
            GenotypeProcessMethod::OrthStandardize, geno);

        REQUIRE(geno.isApproxToConstant(0.0));
    }

    SECTION("all heterozygous")
    {
        // all 1s: maf=0.5
        // OrthCenter recodes: 1→2*0.5=1.0
        // recoded: [1, 1, 1, 1]
        // CenterMethod: mean=1.0, centered: [0, 0, 0, 0]
        // stddev=0, monomorphic, not divided
        Eigen::MatrixXd geno(4, 1);
        geno << 1, 1, 1, 1;

        process_matrix<GeneticEffectType::Dom>(
            GenotypeProcessMethod::OrthStandardize, geno);

        REQUIRE(geno.isApproxToConstant(0.0));
    }
}

// ============================================================================
// OrthStandardizeHWE dominance
// ============================================================================

TEST_CASE("OrthStandardizeHWE dominance mode", "[grm][OrthStandardizeHWE]")
{
    SECTION("basic transformation")
    {
        // genotype: 0, 1, 2, 1, 0
        // maf=0.4, orthogonal Dom recodes: 0→0, 1→0.8, 2→-0.4
        // recoded: [0, 0.8, -0.4, 0.8, 0]
        // OrthHWE mean = 2*maf² = 0.32
        // centered: [-0.32, 0.48, -0.72, 0.48, -0.32]
        // OrthHWE stddev = 2*maf*(1-maf) = 0.48
        // divided by stddev
        Eigen::MatrixXd geno(5, 1);
        geno << 0, 1, 2, 1, 0;

        double stddev = 2.0 * 0.4 * 0.6;

        Eigen::MatrixXd expected(5, 1);
        expected << -0.32 / stddev, 0.48 / stddev, -0.72 / stddev,
            0.48 / stddev, -0.32 / stddev;

        process_matrix<GeneticEffectType::Dom>(
            GenotypeProcessMethod::OrthStandardizeHWE, geno);

        REQUIRE(geno.isApprox(expected));
    }
}

// ============================================================================
// Multiple columns tests
// ============================================================================

TEST_CASE("All policies handle multiple columns", "[grm]")
{
    Eigen::MatrixXd geno(5, 3);
    geno << 0, 1, 2, 1, 2, 0, 2, 1, 1, 1, 0, 2, 0, 2, 1;

    SECTION("Center additive")
    {
        Eigen::MatrixXd g = geno;
        process_matrix<GeneticEffectType::Add>(
            GenotypeProcessMethod::Center, g);
        // verify each column is mean-centered
        for (Eigen::Index i = 0; i < g.cols(); ++i)
        {
            REQUIRE(std::abs(g.col(i).mean()) < 1e-10);
        }
    }

    SECTION("OrthCenter additive")
    {
        Eigen::MatrixXd g = geno;
        process_matrix<GeneticEffectType::Add>(
            GenotypeProcessMethod::OrthCenter, g);
        for (Eigen::Index i = 0; i < g.cols(); ++i)
        {
            REQUIRE(std::abs(g.col(i).mean()) < 1e-10);
        }
    }

    SECTION("OrthStandardize additive - columns standardized")
    {
        Eigen::MatrixXd g = geno;
        process_matrix<GeneticEffectType::Add>(
            GenotypeProcessMethod::OrthStandardize, g);
        // each column should have mean ~ 0
        for (Eigen::Index i = 0; i < g.cols(); ++i)
        {
            REQUIRE(std::abs(g.col(i).mean()) < 1e-10);
        }
    }
}
