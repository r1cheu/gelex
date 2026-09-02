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
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <string>
#include <vector>

#include "gelex/data/bed.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/grm/grm.h"
#include "gelex/data/marker_range.h"
#include "gelex/data/reader.h"
#include "gelex/genetic_mode.h"

#include "bed_fixture.h"

using namespace gelex;  // NOLINT
using Catch::Matchers::WithinAbs;
using gelex::test::are_matrices_equal;
using gelex::test::BedFixture;

// ============================================================================
// GrmBuilder core tests
// ============================================================================

TEST_CASE("GRM - additive GRM", "[data][grm][compute]")
{
    BedFixture fixture;

    SECTION("Happy path - OrthStandardize additive GRM properties")
    {
        const Eigen::Index num_samples = 15;
        const Eigen::Index num_snps = 50;
        auto [bed_prefix, genotypes]
            = fixture.create_bed_files(num_samples, num_snps, 0.0);

        auto bed = open_bed(bed_prefix);
        GrmBuilder builder(
            bed,
            GeneticModeSet{GeneticMode::A},
            GenotypeMethod::OrthStandardize,
            10);
        const std::vector<MarkerRange> ranges{
            {std::string{}, 0, bed.num_snps()}};
        std::vector<GrmMatrix> out;
        builder.build(ranges, [&](const GrmMatrix& m) { out.push_back(m); });
        const GrmMatrix& result = out.at(0);

        // verify dimensions
        REQUIRE(result.grm.rows() == num_samples);
        REQUIRE(result.grm.cols() == num_samples);

        // verify denominator is positive
        REQUIRE(result.denominator > 0.0);

        // Note: cblas_dsyrk only fills lower triangle, upper triangle is 0
        // verify lower triangle has non-zero values (off-diagonal)
        bool has_nonzero_lower = false;
        for (Eigen::Index i = 1; i < num_samples; ++i)
        {
            for (Eigen::Index j = 0; j < i; ++j)
            {
                if (std::abs(result.grm(i, j)) > 1e-10)
                {
                    has_nonzero_lower = true;
                    break;
                }
            }
            if (has_nonzero_lower)
            {
                break;
            }
        }
        REQUIRE(has_nonzero_lower);

        // verify trace normalization: trace / n equals denominator
        double trace_per_n
            = result.grm.trace() / static_cast<double>(num_samples);
        REQUIRE_THAT(trace_per_n, WithinAbs(result.denominator, 1e-10));

        // verify normalized GRM has trace/n = 1.0
        Eigen::MatrixXd G_normalized = result.grm / result.denominator;
        double normalized_trace
            = G_normalized.trace() / static_cast<double>(num_samples);
        REQUIRE_THAT(normalized_trace, WithinAbs(1.0, 1e-10));
    }
}

TEST_CASE("GRM - dominance GRM", "[data][grm][compute]")
{
    BedFixture fixture;

    SECTION("Happy path - OrthStandardize dominance GRM properties")
    {
        const Eigen::Index num_samples = 12;
        const Eigen::Index num_snps = 40;
        auto [bed_prefix, genotypes]
            = fixture.create_bed_files(num_samples, num_snps, 0.0);

        auto bed = open_bed(bed_prefix);
        GrmBuilder builder(
            bed,
            GeneticModeSet{GeneticMode::D},
            GenotypeMethod::OrthStandardize,
            10);
        const std::vector<MarkerRange> ranges{
            {std::string{}, 0, bed.num_snps()}};
        std::vector<GrmMatrix> out;
        builder.build(ranges, [&](const GrmMatrix& m) { out.push_back(m); });
        const GrmMatrix& result = out.at(0);

        // verify dimensions
        REQUIRE(result.grm.rows() == num_samples);
        REQUIRE(result.grm.cols() == num_samples);

        // verify denominator is positive
        REQUIRE(result.denominator > 0.0);

        // Note: cblas_dsyrk only fills lower triangle
        // verify lower triangle has non-zero values
        bool has_nonzero_lower = false;
        for (Eigen::Index i = 1; i < num_samples; ++i)
        {
            for (Eigen::Index j = 0; j < i; ++j)
            {
                if (std::abs(result.grm(i, j)) > 1e-10)
                {
                    has_nonzero_lower = true;
                    break;
                }
            }
            if (has_nonzero_lower)
            {
                break;
            }
        }
        REQUIRE(has_nonzero_lower);

        // verify normalized trace
        Eigen::MatrixXd D_normalized = result.grm / result.denominator;
        double trace_per_n
            = D_normalized.trace() / static_cast<double>(num_samples);
        REQUIRE_THAT(trace_per_n, WithinAbs(1.0, 1e-10));
    }
}

// ============================================================================
// Chunk consistency tests
// ============================================================================

TEST_CASE("GRM - chunk size consistency", "[data][grm][compute][chunk]")
{
    BedFixture fixture;

    SECTION("Different chunk sizes produce identical GRM")
    {
        const Eigen::Index num_samples = 10;
        const Eigen::Index num_snps = 30;
        auto [bed_prefix, genotypes] = fixture.create_bed_files(
            num_samples, num_snps, 0.0, 0.05, 0.5, 42);

        auto bed = open_bed(bed_prefix);
        const auto method = GenotypeMethod::OrthStandardize;
        const std::vector<MarkerRange> ranges{
            {std::string{}, 0, bed.num_snps()}};

        std::vector<GrmMatrix> out;
        auto collect = [&](const GrmMatrix& m) { out.push_back(m); };

        GrmBuilder{bed, GeneticModeSet{GeneticMode::A}, method, 1}.build(
            ranges, collect);
        GrmBuilder{bed, GeneticModeSet{GeneticMode::A}, method, 7}.build(
            ranges, collect);
        GrmBuilder{bed, GeneticModeSet{GeneticMode::A}, method, num_snps}.build(
            ranges, collect);
        GrmBuilder{bed, GeneticModeSet{GeneticMode::A}, method, num_snps + 100}
            .build(ranges, collect);

        REQUIRE(out.size() == 4);
        GrmMatrix& result1 = out.at(0);
        GrmMatrix& result2 = out.at(1);
        GrmMatrix& result3 = out.at(2);
        GrmMatrix& result4 = out.at(3);

        // all should be equal
        REQUIRE(are_matrices_equal(result1.grm, result2.grm, 1e-10));
        REQUIRE(are_matrices_equal(result2.grm, result3.grm, 1e-10));
        REQUIRE(are_matrices_equal(result3.grm, result4.grm, 1e-10));

        // denominators should also be equal
        REQUIRE(std::abs(result1.denominator - result2.denominator) < 1e-10);
        REQUIRE(std::abs(result2.denominator - result3.denominator) < 1e-10);
        REQUIRE(std::abs(result3.denominator - result4.denominator) < 1e-10);
    }
}

// ============================================================================
// Per-chromosome tests
// ============================================================================

TEST_CASE(
    "GRM - per-chromosome GRMs are labelled and split",
    "[data][grm][loco]")
{
    BedFixture fixture;

    SECTION("Two chromosomes yield two labelled GRMs")
    {
        Eigen::MatrixXd genotypes{
            {0, 1, 2, 1}, {1, 1, 1, 0}, {2, 1, 0, 2}, {1, 0, 1, 1}};

        // two markers on chr1, two on chr2
        auto [bed_prefix, _] = fixture.create_deterministic_bed_files(
            genotypes, {}, {}, {"1", "1", "2", "2"});

        auto bed = open_bed(bed_prefix);
        auto bim = read_bim(bed_prefix.string() + ".bim");
        auto ranges = chromosome_ranges(bim);
        REQUIRE(ranges.size() == 2);

        GrmBuilder builder(
            bed, GeneticModeSet{GeneticMode::A}, GenotypeMethod::Center, 10);
        std::vector<GrmMatrix> out;
        builder.build(ranges, [&](const GrmMatrix& m) { out.push_back(m); });

        REQUIRE(out.size() == 2);
        REQUIRE(out.at(0).label == "1");
        REQUIRE(out.at(1).label == "2");
        REQUIRE(out.at(0).mode == GeneticMode::A);
    }
}

// ============================================================================
// Numerical correctness tests
// ============================================================================

TEST_CASE("GRM - numerical correctness", "[data][grm][compute][numerical]")
{
    BedFixture fixture;

    SECTION("OrthStandardize additive GRM with deterministic genotype")
    {
        // create a simple deterministic genotype matrix
        // 4 samples, 3 SNPs
        Eigen::MatrixXd genotypes{{0, 1, 2}, {1, 1, 1}, {2, 1, 0}, {1, 1, 1}};

        auto [bed_prefix, _]
            = fixture.create_deterministic_bed_files(genotypes);

        auto bed = open_bed(bed_prefix);
        GrmBuilder builder(
            bed,
            GeneticModeSet{GeneticMode::A},
            GenotypeMethod::OrthStandardize,
            10);
        const std::vector<MarkerRange> ranges{
            {std::string{}, 0, bed.num_snps()}};
        std::vector<GrmMatrix> out;
        builder.build(ranges, [&](const GrmMatrix& m) { out.push_back(m); });
        const GrmMatrix& result = out.at(0);

        // manually compute expected GRM using OrthStandardize additive method
        // OrthStandardizeMethod = CenterMethod + divide by sample stddev
        Eigen::MatrixXd Z = genotypes;
        for (Eigen::Index j = 0; j < Z.cols(); ++j)
        {
            double mean = Z.col(j).mean();
            Z.col(j).array() -= mean;
            double var = Z.col(j).squaredNorm() / static_cast<double>(Z.rows());
            double denom = std::sqrt(var);
            if (denom > 1e-10)
            {
                Z.col(j).array() /= denom;
            }
            else
            {
                Z.col(j).setZero();
            }
        }

        // Step 2: compute GRM = Z * Z^T (unnormalized)
        Eigen::MatrixXd expected_unnormalized_grm = Z * Z.transpose();

        // Step 3: compute expected denominator
        auto n = static_cast<double>(genotypes.rows());
        double expected_denom = expected_unnormalized_grm.trace() / n;

        // Verify denominator matches
        REQUIRE_THAT(result.denominator, WithinAbs(expected_denom, 1e-8));

        // Step 4: normalize for comparison
        Eigen::MatrixXd expected_grm
            = expected_unnormalized_grm / expected_denom;
        Eigen::MatrixXd G = result.grm / result.denominator;

        // verify diagonal elements
        for (Eigen::Index i = 0; i < 4; ++i)
        {
            REQUIRE_THAT(G(i, i), WithinAbs(expected_grm(i, i), 1e-8));
        }

        // verify lower triangle elements (dsyrk only fills lower triangle)
        for (Eigen::Index i = 1; i < 4; ++i)
        {
            for (Eigen::Index j = 0; j < i; ++j)
            {
                REQUIRE_THAT(G(i, j), WithinAbs(expected_grm(i, j), 1e-8));
            }
        }
    }

    SECTION("Center additive GRM with deterministic genotype")
    {
        // 3 samples, 2 SNPs
        Eigen::MatrixXd genotypes{{0, 2}, {1, 1}, {2, 0}};

        auto [bed_prefix, _]
            = fixture.create_deterministic_bed_files(genotypes);

        auto bed = open_bed(bed_prefix);
        GrmBuilder builder(
            bed, GeneticModeSet{GeneticMode::A}, GenotypeMethod::Center, 10);
        const std::vector<MarkerRange> ranges{
            {std::string{}, 0, bed.num_snps()}};
        std::vector<GrmMatrix> out;
        builder.build(ranges, [&](const GrmMatrix& m) { out.push_back(m); });
        const GrmMatrix& result = out.at(0);

        // Center additive: mean centering
        Eigen::MatrixXd Z = genotypes;
        for (Eigen::Index j = 0; j < Z.cols(); ++j)
        {
            Z.col(j).array() -= Z.col(j).mean();
        }

        Eigen::MatrixXd expected_unnormalized = Z * Z.transpose();
        auto n = static_cast<double>(genotypes.rows());
        double expected_denom = expected_unnormalized.trace() / n;

        // Verify denominator
        REQUIRE_THAT(result.denominator, WithinAbs(expected_denom, 1e-8));

        // Normalize for comparison
        Eigen::MatrixXd expected_grm = expected_unnormalized / expected_denom;
        Eigen::MatrixXd G = result.grm / result.denominator;

        // verify diagonal elements
        for (Eigen::Index i = 0; i < 3; ++i)
        {
            REQUIRE_THAT(G(i, i), WithinAbs(expected_grm(i, i), 1e-8));
        }

        // verify lower triangle
        for (Eigen::Index i = 1; i < 3; ++i)
        {
            for (Eigen::Index j = 0; j < i; ++j)
            {
                REQUIRE_THAT(G(i, j), WithinAbs(expected_grm(i, j), 1e-8));
            }
        }
    }
}
