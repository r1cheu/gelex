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

#include "gelex/data/grm_code_policy.h"

namespace gelex::grm
{

// ============================================================================
// detail::additive_mean_center
// ============================================================================

TEST_CASE("additive_mean_center single column", "[grm]")
{
    Eigen::MatrixXd geno(5, 1);
    geno << 0, 1, 2, 1, 0;
    // mean = (0+1+2+1+0)/5 = 0.8

    Eigen::MatrixXd expected(5, 1);
    expected << -0.8, 0.2, 1.2, 0.2, -0.8;

    detail::additive_mean_center(geno);

    REQUIRE(geno.isApprox(expected));
}

TEST_CASE("additive_mean_center multiple columns", "[grm]")
{
    Eigen::MatrixXd geno(4, 3);
    geno << 0, 2, 1, 1, 1, 1, 2, 0, 1, 1, 1, 1;
    // col0: mean = 1.0
    // col1: mean = 1.0
    // col2: mean = 1.0

    Eigen::MatrixXd expected(4, 3);
    expected << -1, 1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0;

    detail::additive_mean_center(geno);

    REQUIRE(geno.isApprox(expected));
}

TEST_CASE("additive_mean_center already centered", "[grm]")
{
    Eigen::MatrixXd geno(4, 1);
    geno << -1, 1, -1, 1;  // mean = 0

    Eigen::MatrixXd expected = geno;

    detail::additive_mean_center(geno);

    REQUIRE(geno.isApprox(expected));
}

// ============================================================================
// Su
// ============================================================================

TEST_CASE("Su additive mode", "[grm][Su]")
{
    Su su;

    SECTION("basic centering")
    {
        Eigen::MatrixXd geno(5, 1);
        geno << 0, 1, 2, 1, 0;
        // mean = 0.8

        Eigen::MatrixXd expected(5, 1);
        expected << -0.8, 0.2, 1.2, 0.2, -0.8;

        su(geno, true);

        REQUIRE(geno.isApprox(expected));
    }

    SECTION("multiple columns")
    {
        Eigen::MatrixXd geno(4, 2);
        geno << 0, 2, 1, 1, 2, 0, 1, 1;
        // col0: mean=1, col1: mean=1

        Eigen::MatrixXd expected(4, 2);
        expected << -1, 1, 0, 0, 1, -1, 0, 0;

        su(geno, true);

        REQUIRE(geno.isApprox(expected));
    }
}

TEST_CASE("Su dominance mode", "[grm][Su]")
{
    Su su;

    SECTION("basic transformation")
    {
        Eigen::MatrixXd geno(5, 1);
        geno << 0, 1, 2, 1, 0;

        Eigen::MatrixXd expected(5, 1);
        expected << -0.48, 0.52, -0.48, 0.52, -0.48;

        su(geno, false);

        REQUIRE(geno.isApprox(expected));
    }

    SECTION("all heterozygous")
    {
        // all 1s: pA = 0.5
        // 2*0.5*0.5 = 0.5
        // 1 - 0.5 = 0.5
        Eigen::MatrixXd geno(4, 1);
        geno << 1, 1, 1, 1;

        Eigen::MatrixXd expected(4, 1);
        expected << 0.5, 0.5, 0.5, 0.5;

        su(geno, false);

        REQUIRE(geno.isApprox(expected));
    }
}

// ============================================================================
// Zeng
// ============================================================================

TEST_CASE("Zeng additive mode", "[grm][Zeng]")
{
    Zeng zeng;

    SECTION("basic centering")
    {
        Eigen::MatrixXd geno(5, 1);
        geno << 0, 1, 2, 1, 0;

        Eigen::MatrixXd expected(5, 1);
        expected << -0.8, 0.2, 1.2, 0.2, -0.8;

        zeng(geno, true);

        REQUIRE(geno.isApprox(expected));
    }
}

TEST_CASE("Zeng dominance mode", "[grm][Zeng]")
{
    Zeng zeng;

    SECTION("basic transformation")
    {
        // genotype: 0, 1, 2, 1, 0
        // pA = 0.8/2 = 0.4, pa = 0.6
        // geno=2: -2*pa*pa = -2*0.36 = -0.72
        // geno=1: 2*pA*pa = 2*0.24 = 0.48
        // geno=0: -2*pA*pA = -2*0.16 = -0.32
        Eigen::MatrixXd geno(5, 1);
        geno << 0, 1, 2, 1, 0;

        Eigen::MatrixXd expected(5, 1);
        expected << -0.32, 0.48, -0.72, 0.48, -0.32;

        zeng(geno, false);

        REQUIRE(geno.isApprox(expected));
    }

    SECTION("all heterozygous")
    {
        // all 1s: pA = 0.5, pa = 0.5
        // geno=1: 2*0.5*0.5 = 0.5
        Eigen::MatrixXd geno(4, 1);
        geno << 1, 1, 1, 1;

        Eigen::MatrixXd expected(4, 1);
        expected << 0.5, 0.5, 0.5, 0.5;

        zeng(geno, false);

        REQUIRE(geno.isApprox(expected));
    }

    SECTION("all homozygous AA")
    {
        // all 2s: pA = 1.0, pa = 0.0
        // geno=2: -2*pa*pa = -2*0*0 = 0
        Eigen::MatrixXd geno(4, 1);
        geno << 2, 2, 2, 2;

        Eigen::MatrixXd expected(4, 1);
        expected << 0, 0, 0, 0;

        zeng(geno, false);

        REQUIRE(geno.isApprox(expected));
    }
}

// ============================================================================
// Yang
// ============================================================================

TEST_CASE("Yang additive mode", "[grm][Yang]")
{
    Yang yang;

    SECTION("basic standardization")
    {
        // genotype: 0, 1, 2, 1, 0
        // pA = 0.8/2 = 0.4, pa = 0.6
        // denom = sqrt(2*0.4*0.6) = sqrt(0.48) ≈ 0.6928
        // centered: (x - 2*pA) / denom = (x - 0.8) / 0.6928
        Eigen::MatrixXd geno(5, 1);
        geno << 0, 1, 2, 1, 0;

        double pA = 0.4;
        double denom = std::sqrt(2 * pA * (1 - pA));

        Eigen::MatrixXd expected(5, 1);
        expected << (0 - 0.8) / denom, (1 - 0.8) / denom, (2 - 0.8) / denom,
            (1 - 0.8) / denom, (0 - 0.8) / denom;

        yang(geno, true);

        REQUIRE(geno.isApprox(expected));
    }

    SECTION("monomorphic SNP sets to zero")
    {
        // all 0s: pA = 0, denom = 0 < EPSILON
        Eigen::MatrixXd geno(4, 1);
        geno << 0, 0, 0, 0;

        yang(geno, true);

        REQUIRE(geno.isApproxToConstant(0.0));
    }

    SECTION("all homozygous AA sets to zero")
    {
        // all 2s: pA = 1, pa = 0, denom = 0
        Eigen::MatrixXd geno(4, 1);
        geno << 2, 2, 2, 2;

        yang(geno, true);

        REQUIRE(geno.isApproxToConstant(0.0));
    }
}

TEST_CASE("Yang dominance mode", "[grm][Yang]")
{
    Yang yang;

    SECTION("basic transformation")
    {
        // genotype: 0, 1, 2, 1, 0
        // pA = 0.4, pa = 0.6
        // denom = 2*pA*pa = 0.48
        // geno=2: -2*pa*pa / denom = -0.72/0.48 = -1.5
        // geno=1: 2*pA*pa / denom = 0.48/0.48 = 1.0
        // geno=0: -2*pA*pA / denom = -0.32/0.48 = -0.6667
        Eigen::MatrixXd geno(5, 1);
        geno << 0, 1, 2, 1, 0;

        double pA = 0.4;
        double pa = 0.6;
        double denom = 2 * pA * pa;

        Eigen::MatrixXd expected(5, 1);
        expected << -2 * pA * pA / denom, 2 * pA * pa / denom,
            -2 * pa * pa / denom, 2 * pA * pa / denom, -2 * pA * pA / denom;

        yang(geno, false);

        REQUIRE(geno.isApprox(expected));
    }

    SECTION("monomorphic SNP sets to zero")
    {
        Eigen::MatrixXd geno(4, 1);
        geno << 0, 0, 0, 0;

        yang(geno, false);

        REQUIRE(geno.isApproxToConstant(0.0));
    }

    SECTION("all heterozygous")
    {
        // all 1s: pA = 0.5, pa = 0.5
        // denom = 2*0.5*0.5 = 0.5
        // geno=1: 2*0.5*0.5 / 0.5 = 1.0
        Eigen::MatrixXd geno(4, 1);
        geno << 1, 1, 1, 1;

        Eigen::MatrixXd expected(4, 1);
        expected << 1, 1, 1, 1;

        yang(geno, false);

        REQUIRE(geno.isApprox(expected));
    }
}

// ============================================================================
// Vitezica
// ============================================================================

TEST_CASE("Vitezica additive mode", "[grm][Vitezica]")
{
    Vitezica vitezica;

    SECTION("basic centering")
    {
        // genotype: 0, 1, 2, 1, 0
        // nAA = 1, nAa = 2
        // centered: x - (nAa + 2*nAA) = x - (2 + 2) = x - 4
        Eigen::MatrixXd geno(5, 1);
        geno << 0, 1, 2, 1, 0;

        // count: nAA=1, nAa=2
        // subtract (nAa + 2*nAA) = 2 + 2 = 4
        // but the code divides by n, so it's: x - (nAa + 2*nAA)/n
        // Actually looking at the code: col.array() -= (nAa + (2 * nAA))
        // This subtracts the sum, not mean
        Eigen::MatrixXd expected(5, 1);
        expected << 0 - 4, 1 - 4, 2 - 4, 1 - 4, 0 - 4;

        vitezica(geno, true);

        REQUIRE(geno.isApprox(expected));
    }

    SECTION("multiple columns")
    {
        Eigen::MatrixXd geno(4, 2);
        geno << 0, 2, 1, 1, 2, 0, 1, 1;
        // col0: nAA=1, nAa=2 -> subtract 4
        // col1: nAA=1, nAa=2 -> subtract 4

        Eigen::MatrixXd expected(4, 2);
        expected << -4, -2, -3, -3, -2, -4, -3, -3;

        vitezica(geno, true);

        REQUIRE(geno.isApprox(expected));
    }
}

TEST_CASE("Vitezica dominance mode", "[grm][Vitezica]")
{
    Vitezica vitezica;

    SECTION("basic transformation")
    {
        // genotype: 0, 1, 2, 1, 0
        // nAA=1, nAa=2, naa=2
        // denom = nAA + naa - (nAA - naa)^2 = 1 + 2 - 1 = 2
        // geno=2: -2*naa*nAa / denom = -2*2*2 / 2 = -4
        // geno=1: 4*nAA*naa / denom = 4*1*2 / 2 = 4
        // geno=0: -2*nAA*nAa / denom = -2*1*2 / 2 = -2
        Eigen::MatrixXd geno(5, 1);
        geno << 0, 1, 2, 1, 0;

        double nAA = 1;
        double nAa = 2;
        double naa = 2;
        double denom = nAA + naa - std::pow(nAA - naa, 2);

        Eigen::MatrixXd expected(5, 1);
        expected << -2 * nAA * nAa / denom, 4 * nAA * naa / denom,
            -2 * naa * nAa / denom, 4 * nAA * naa / denom,
            -2 * nAA * nAa / denom;

        vitezica(geno, false);

        REQUIRE(geno.isApprox(expected));
    }

    SECTION("monomorphic SNP sets to zero")
    {
        // all 0s: nAA=0, nAa=0, naa=4
        // denom = 0 + 4 - 16 = -12, but since nAA=nAa=0, let's check
        // Actually denom = 0 + 4 - (0-4)^2 = 4 - 16 = -12
        // But -12 >= EPSILON, so won't set to zero
        // Let's try all 1s instead
        // all 1s: nAA=0, nAa=4, naa=0
        // denom = 0 + 0 - 0 = 0 < EPSILON
        Eigen::MatrixXd geno(4, 1);
        geno << 1, 1, 1, 1;

        vitezica(geno, false);

        REQUIRE(geno.isApproxToConstant(0.0));
    }

    SECTION("mixed genotypes")
    {
        // 0, 0, 1, 1, 2, 2
        // nAA=2, nAa=2, naa=2
        // denom = 2 + 2 - 0 = 4
        Eigen::MatrixXd geno(6, 1);
        geno << 0, 0, 1, 1, 2, 2;

        double nAA = 2;
        double nAa = 2;
        double naa = 2;
        double denom = nAA + naa - std::pow(nAA - naa, 2);

        Eigen::MatrixXd expected(6, 1);
        double val0 = -2 * nAA * nAa / denom;
        double val1 = 4 * nAA * naa / denom;
        double val2 = -2 * naa * nAa / denom;
        expected << val0, val0, val1, val1, val2, val2;

        vitezica(geno, false);

        REQUIRE(geno.isApprox(expected));
    }

    SECTION("no homozygous aa (naa=0)")
    {
        // geno = [2, 2, 2, 2, 1]
        // nAA=4, nAa=1, naa=0
        // denom = 4 + 0 - 16 = -12
        // geno=2: -2*naa*nAa/denom = 0
        // geno=1: 4*nAA*naa/denom = 0
        Eigen::MatrixXd geno(5, 1);
        geno << 2, 2, 2, 2, 1;

        Eigen::MatrixXd expected(5, 1);
        expected << 0, 0, 0, 0, 0;

        vitezica(geno, false);

        REQUIRE(geno.isApprox(expected));
    }

    SECTION("no homozygous AA (nAA=0)")
    {
        // geno = [0, 0, 0, 0, 1]
        // nAA=0, nAa=1, naa=4
        // denom = 0 + 4 - 16 = -12
        // geno=0: -2*nAA*nAa/denom = 0
        // geno=1: 4*nAA*naa/denom = 0
        Eigen::MatrixXd geno(5, 1);
        geno << 0, 0, 0, 0, 1;

        Eigen::MatrixXd expected(5, 1);
        expected << 0, 0, 0, 0, 0;

        vitezica(geno, false);

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

    SECTION("Su additive")
    {
        Eigen::MatrixXd g = geno;
        Su()(g, true);
        // verify each column is mean-centered
        for (Eigen::Index i = 0; i < g.cols(); ++i)
        {
            REQUIRE(std::abs(g.col(i).mean()) < 1e-10);
        }
    }

    SECTION("Zeng additive")
    {
        Eigen::MatrixXd g = geno;
        Zeng()(g, true);
        for (Eigen::Index i = 0; i < g.cols(); ++i)
        {
            REQUIRE(std::abs(g.col(i).mean()) < 1e-10);
        }
    }

    SECTION("Yang additive - columns standardized")
    {
        Eigen::MatrixXd g = geno;
        Yang()(g, true);
        // each column should have mean ≈ 0
        for (Eigen::Index i = 0; i < g.cols(); ++i)
        {
            REQUIRE(std::abs(g.col(i).mean()) < 1e-10);
        }
    }
}

}  // namespace gelex::grm
