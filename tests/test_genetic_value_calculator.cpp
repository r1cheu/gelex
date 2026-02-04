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
#include <unordered_map>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "gelex/data/effect_sampler.h"
#include "gelex/data/genetic_value_calculator.h"

using namespace gelex;  // NOLINT
using Catch::Matchers::WithinAbs;

TEST_CASE("GeneticValueCalculator - basic calculation", "[genetic_calc]")
{
    SECTION("Single SNP with unit effect")
    {
        // 3 individuals, 1 SNP
        Eigen::MatrixXd geno(3, 1);
        geno << 0, 1, 2;

        std::unordered_map<Eigen::Index, CausalEffect> effects{
            {0, {.additive = 1.0, .dominance = 0.0, .add_class = 0, .dom_class = 0}},
        };

        auto result = GeneticValueCalculator::calculate_chunk(
            geno, effects, 0, geno.cols(), false);

        REQUIRE(result.additive.size() == 3);
        REQUIRE(result.dominance.size() == 3);

        // Dominance values should be zero
        for (Eigen::Index i = 0; i < 3; ++i)
        {
            REQUIRE(result.dominance(i) == 0.0);
        }
    }

    SECTION("Multiple SNPs with different effects")
    {
        // 4 individuals, 3 SNPs
        Eigen::MatrixXd geno(4, 3);
        geno << 0, 1, 2, 1, 1, 1, 2, 0, 0, 1, 2, 1;

        std::unordered_map<Eigen::Index, CausalEffect> effects{
            {0, {.additive = 1.0, .dominance = 0.0, .add_class = 0, .dom_class = 0}},
            {1, {.additive = 0.5, .dominance = 0.0, .add_class = 0, .dom_class = 0}},
            {2, {.additive = -0.5, .dominance = 0.0, .add_class = 0, .dom_class = 0}},
        };

        auto result = GeneticValueCalculator::calculate_chunk(
            geno, effects, 0, geno.cols(), false);

        REQUIRE(result.additive.size() == 4);
    }

    SECTION("Only subset of SNPs have effects")
    {
        Eigen::MatrixXd geno(3, 5);
        geno.setOnes();

        // Only SNP 1 and SNP 3 have effects
        std::unordered_map<Eigen::Index, CausalEffect> effects{
            {1, {.additive = 1.0, .dominance = 0.0, .add_class = 0, .dom_class = 0}},
            {3, {.additive = 2.0, .dominance = 0.0, .add_class = 0, .dom_class = 0}},
        };

        auto result = GeneticValueCalculator::calculate_chunk(
            geno, effects, 0, geno.cols(), false);

        REQUIRE(result.additive.size() == 3);
    }

    SECTION("Empty effects produces zero genetic values")
    {
        Eigen::MatrixXd geno(3, 5);
        geno.setRandom();

        std::unordered_map<Eigen::Index, CausalEffect> effects;

        auto result = GeneticValueCalculator::calculate_chunk(
            geno, effects, 0, geno.cols(), false);

        for (Eigen::Index i = 0; i < 3; ++i)
        {
            REQUIRE(result.additive(i) == 0.0);
            REQUIRE(result.dominance(i) == 0.0);
        }
    }
}

TEST_CASE("GeneticValueCalculator - dominance effects", "[genetic_calc]")
{
    SECTION("Dominance values computed when has_dominance is true")
    {
        Eigen::MatrixXd geno(3, 2);
        geno << 0, 2, 1, 1, 2, 0;

        std::unordered_map<Eigen::Index, CausalEffect> effects{
            {0,
             {.additive = 1.0, .dominance = 0.5, .add_class = 0, .dom_class = 0}},
            {1,
             {.additive = 0.5, .dominance = 1.0, .add_class = 0, .dom_class = 0}},
        };

        auto result = GeneticValueCalculator::calculate_chunk(
            geno, effects, 0, geno.cols(), true);

        REQUIRE(result.additive.size() == 3);
        REQUIRE(result.dominance.size() == 3);

        // Check that dominance values are not all zero
        bool has_nonzero = false;
        for (Eigen::Index i = 0; i < 3; ++i)
        {
            if (std::abs(result.dominance(i)) > 1e-10)
            {
                has_nonzero = true;
                break;
            }
        }
        REQUIRE(has_nonzero);
    }

    SECTION("Dominance values are zero when has_dominance is false")
    {
        Eigen::MatrixXd geno(3, 2);
        geno << 0, 2, 1, 1, 2, 0;

        std::unordered_map<Eigen::Index, CausalEffect> effects{
            {0,
             {.additive = 1.0, .dominance = 0.5, .add_class = 0, .dom_class = 0}},
        };

        auto result = GeneticValueCalculator::calculate_chunk(
            geno, effects, 0, geno.cols(), false);

        for (Eigen::Index i = 0; i < 3; ++i)
        {
            REQUIRE(result.dominance(i) == 0.0);
        }
    }
}

TEST_CASE("GeneticValueCalculator - chunk calculation", "[genetic_calc]")
{
    SECTION("Chunk calculation filters by index range")
    {
        // Use varying genotypes to avoid all-zero after standardization
        Eigen::MatrixXd chunk(3, 5);
        chunk << 0, 1, 2, 1, 0, 1, 1, 1, 2, 1, 2, 0, 0, 0, 2;

        // Global indices 10-14, only 12 has effect
        std::unordered_map<Eigen::Index, CausalEffect> effects{
            {5,
             {.additive = 100.0,
              .dominance = 0.0,
              .add_class = 0,
              .dom_class = 0}},  // outside range
            {12,
             {.additive = 1.0,
              .dominance = 0.0,
              .add_class = 0,
              .dom_class = 0}},  // inside range (col 2)
            {20,
             {.additive = 100.0,
              .dominance = 0.0,
              .add_class = 0,
              .dom_class = 0}},  // outside range
        };

        auto result = GeneticValueCalculator::calculate_chunk(
            chunk, effects, 10, 15, false);

        REQUIRE(result.additive.size() == 3);

        // Values should be non-zero only due to effect at index 12
        // Not all zeros (because SNP 12 has effect with varying genotypes)
        bool all_zero = true;
        for (Eigen::Index i = 0; i < 3; ++i)
        {
            if (std::abs(result.additive(i)) > 1e-10)
            {
                all_zero = false;
                break;
            }
        }
        REQUIRE_FALSE(all_zero);
    }

    SECTION("No effects in range produces zero values")
    {
        Eigen::MatrixXd chunk(3, 5);
        chunk.setRandom();

        std::unordered_map<Eigen::Index, CausalEffect> effects{
            {0, {.additive = 1.0, .dominance = 0.0, .add_class = 0, .dom_class = 0}},
            {100, {.additive = 1.0, .dominance = 0.0, .add_class = 0, .dom_class = 0}},
        };

        auto result = GeneticValueCalculator::calculate_chunk(
            chunk, effects, 10, 15, false);

        for (Eigen::Index i = 0; i < 3; ++i)
        {
            REQUIRE(result.additive(i) == 0.0);
        }
    }
}
