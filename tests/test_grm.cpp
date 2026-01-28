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
#include <filesystem>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "bed_fixture.h"
#include "gelex/data/grm.h"
#include "gelex/data/grm_code_policy.h"

namespace fs = std::filesystem;

using namespace gelex;  // NOLINT
using Catch::Matchers::WithinAbs;
using gelex::test::are_matrices_equal;
using gelex::test::BedFixture;

// ============================================================================
// Construction tests
// ============================================================================

TEST_CASE("GRM - Construction with valid BED files", "[grm][construction]")
{
    BedFixture fixture;

    SECTION("Happy path - construct from valid BED prefix")
    {
        const Eigen::Index num_samples = 10;
        const Eigen::Index num_snps = 20;
        auto [bed_prefix, genotypes]
            = fixture.create_bed_files(num_samples, num_snps, 0.0);

        REQUIRE_NOTHROW(
            [&]()
            {
                GRM grm(bed_prefix);
                REQUIRE(grm.num_snps() == num_snps);
                REQUIRE(
                    grm.sample_ids().size()
                    == static_cast<size_t>(num_samples));
            }());
    }

    SECTION("Exception - file not found")
    {
        fs::path non_existent_path = "/tmp/non_existent_file_12345";

        // GRM constructor first creates SampleManager from FAM file,
        // so exception is thrown for missing FAM file
        REQUIRE_THROWS(GRM(non_existent_path));
    }
}

// ============================================================================
// compute() core tests
// ============================================================================

TEST_CASE("GRM - compute() additive GRM", "[grm][compute]")
{
    BedFixture fixture;

    SECTION("Happy path - Yang additive GRM properties")
    {
        const Eigen::Index num_samples = 15;
        const Eigen::Index num_snps = 50;
        auto [bed_prefix, genotypes]
            = fixture.create_bed_files(num_samples, num_snps, 0.0);

        GRM grm(bed_prefix);
        GrmResult result
            = grm.compute<grm::Yang>(10, true);  // additive, chunk_size=10

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

TEST_CASE("GRM - compute() dominance GRM", "[grm][compute]")
{
    BedFixture fixture;

    SECTION("Happy path - Yang dominance GRM properties")
    {
        const Eigen::Index num_samples = 12;
        const Eigen::Index num_snps = 40;
        auto [bed_prefix, genotypes]
            = fixture.create_bed_files(num_samples, num_snps, 0.0);

        GRM grm(bed_prefix);
        GrmResult result
            = grm.compute<grm::Yang>(10, false);  // dominance, chunk_size=10

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

TEST_CASE("GRM - compute() chunk size consistency", "[grm][compute][chunk]")
{
    BedFixture fixture;

    SECTION("Different chunk sizes produce identical GRM")
    {
        const Eigen::Index num_samples = 10;
        const Eigen::Index num_snps = 30;
        auto [bed_prefix, genotypes] = fixture.create_bed_files(
            num_samples, num_snps, 0.0, 0.05, 0.5, 42);

        // compute with different chunk sizes using same data
        GRM grm1(bed_prefix);
        GrmResult result1
            = grm1.compute<grm::Yang>(1, true);  // chunk_size=1

        GRM grm2(bed_prefix);
        GrmResult result2
            = grm2.compute<grm::Yang>(7, true);  // chunk_size=7

        GRM grm3(bed_prefix);
        GrmResult result3
            = grm3.compute<grm::Yang>(num_snps, true);  // chunk_size=all

        GRM grm4(bed_prefix);
        GrmResult result4
            = grm4.compute<grm::Yang>(num_snps + 100, true);  // chunk_size > n

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
// Accessor tests
// ============================================================================

TEST_CASE("GRM - accessor methods", "[grm][accessor]")
{
    BedFixture fixture;

    SECTION("sample_ids() returns correct sample IDs")
    {
        const Eigen::Index num_samples = 5;
        const Eigen::Index num_snps = 10;
        auto [bed_prefix, genotypes]
            = fixture.create_bed_files(num_samples, num_snps, 0.0);

        GRM grm(bed_prefix);

        const auto& ids = grm.sample_ids();
        REQUIRE(ids.size() == static_cast<size_t>(num_samples));

        // IDs should not be empty
        for (const auto& id : ids)
        {
            REQUIRE_FALSE(id.empty());
        }
    }

    SECTION("num_snps() returns correct count")
    {
        const Eigen::Index num_samples = 8;
        const Eigen::Index num_snps = 25;
        auto [bed_prefix, genotypes]
            = fixture.create_bed_files(num_samples, num_snps, 0.0);

        GRM grm(bed_prefix);

        REQUIRE(grm.num_snps() == num_snps);
    }
}

// ============================================================================
// Numerical correctness tests
// ============================================================================

TEST_CASE("GRM - numerical correctness", "[grm][compute][numerical]")
{
    BedFixture fixture;

    SECTION("Yang additive GRM with deterministic genotype")
    {
        // create a simple deterministic genotype matrix
        // 4 samples, 3 SNPs
        Eigen::MatrixXd genotypes(4, 3);
        // clang-format off
        genotypes << 0, 1, 2,
                     1, 1, 1,
                     2, 1, 0,
                     1, 1, 1;
        // clang-format on

        auto [bed_prefix, _]
            = fixture.create_deterministic_bed_files(genotypes);

        GRM grm(bed_prefix);
        GrmResult result = grm.compute<grm::Yang>(10, true);

        // manually compute expected GRM using Yang additive method
        // Step 1: standardize each column
        Eigen::MatrixXd Z = genotypes;
        for (Eigen::Index j = 0; j < Z.cols(); ++j)
        {
            double pA = Z.col(j).mean() / 2.0;
            double pa = 1.0 - pA;
            double denom = std::sqrt(2.0 * pA * pa);
            if (denom > 1e-10)
            {
                Z.col(j).array() -= 2.0 * pA;
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
        Eigen::MatrixXd expected_grm = expected_unnormalized_grm / expected_denom;
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

    SECTION("Su additive GRM (mean-centered) with deterministic genotype")
    {
        // 3 samples, 2 SNPs
        Eigen::MatrixXd genotypes(3, 2);
        // clang-format off
        genotypes << 0, 2,
                     1, 1,
                     2, 0;
        // clang-format on

        auto [bed_prefix, _]
            = fixture.create_deterministic_bed_files(genotypes);

        GRM grm(bed_prefix);
        GrmResult result = grm.compute<grm::Su>(10, true);

        // Su additive: mean centering
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

// ============================================================================
// Progress callback tests
// ============================================================================

TEST_CASE("GRM - progress callback", "[grm][compute][callback]")
{
    BedFixture fixture;

    SECTION("Progress callback is invoked correctly")
    {
        const Eigen::Index num_samples = 6;
        const Eigen::Index num_snps = 20;
        auto [bed_prefix, genotypes]
            = fixture.create_bed_files(num_samples, num_snps, 0.0);

        GRM grm(bed_prefix);

        std::vector<std::pair<Eigen::Index, Eigen::Index>> progress_calls;

        auto callback
            = [&progress_calls](Eigen::Index current, Eigen::Index total)
        { progress_calls.emplace_back(current, total); };

        GrmResult result = grm.compute<grm::Yang>(5, true, callback);

        // with chunk_size=5 and num_snps=20, expect 4 callback calls
        REQUIRE(progress_calls.size() == 4);

        // verify total is always num_snps
        for (const auto& [current, total] : progress_calls)
        {
            REQUIRE(total == num_snps);
        }

        // verify progress increases
        REQUIRE(progress_calls[0].first == 5);
        REQUIRE(progress_calls[1].first == 10);
        REQUIRE(progress_calls[2].first == 15);
        REQUIRE(progress_calls[3].first == 20);

        // verify result is valid
        REQUIRE(result.grm.rows() == num_samples);
        REQUIRE(result.denominator > 0.0);
    }

    SECTION("Null callback does not crash")
    {
        const Eigen::Index num_samples = 4;
        const Eigen::Index num_snps = 10;
        auto [bed_prefix, genotypes]
            = fixture.create_bed_files(num_samples, num_snps, 0.0);

        GRM grm(bed_prefix);

        GrmResult result = grm.compute<grm::Yang>(5, true, nullptr);
        REQUIRE(result.grm.rows() == num_samples);
        REQUIRE(result.denominator > 0.0);
    }
}
