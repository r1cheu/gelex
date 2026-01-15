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
        Eigen::MatrixXd G
            = grm.compute<grm::Yang>(10, true);  // additive, chunk_size=10

        // verify dimensions
        REQUIRE(G.rows() == num_samples);
        REQUIRE(G.cols() == num_samples);

        // Note: cblas_dsyrk only fills lower triangle, upper triangle is 0
        // verify lower triangle has non-zero values (off-diagonal)
        bool has_nonzero_lower = false;
        for (Eigen::Index i = 1; i < num_samples; ++i)
        {
            for (Eigen::Index j = 0; j < i; ++j)
            {
                if (std::abs(G(i, j)) > 1e-10)
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

        // verify trace normalization: trace / n = 1.0
        double trace_per_n = G.trace() / static_cast<double>(num_samples);
        REQUIRE_THAT(trace_per_n, WithinAbs(1.0, 1e-10));
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
        Eigen::MatrixXd D
            = grm.compute<grm::Yang>(10, false);  // dominance, chunk_size=10

        // verify dimensions
        REQUIRE(D.rows() == num_samples);
        REQUIRE(D.cols() == num_samples);

        // Note: cblas_dsyrk only fills lower triangle
        // verify lower triangle has non-zero values
        bool has_nonzero_lower = false;
        for (Eigen::Index i = 1; i < num_samples; ++i)
        {
            for (Eigen::Index j = 0; j < i; ++j)
            {
                if (std::abs(D(i, j)) > 1e-10)
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

        // verify trace normalization
        double trace_per_n = D.trace() / static_cast<double>(num_samples);
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
        Eigen::MatrixXd G1 = grm1.compute<grm::Yang>(1, true);  // chunk_size=1

        GRM grm2(bed_prefix);
        Eigen::MatrixXd G2 = grm2.compute<grm::Yang>(7, true);  // chunk_size=7

        GRM grm3(bed_prefix);
        Eigen::MatrixXd G3
            = grm3.compute<grm::Yang>(num_snps, true);  // chunk_size=all

        GRM grm4(bed_prefix);
        Eigen::MatrixXd G4
            = grm4.compute<grm::Yang>(num_snps + 100, true);  // chunk_size > n

        // all should be equal
        REQUIRE(are_matrices_equal(G1, G2, 1e-10));
        REQUIRE(are_matrices_equal(G2, G3, 1e-10));
        REQUIRE(are_matrices_equal(G3, G4, 1e-10));
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
        Eigen::MatrixXd G = grm.compute<grm::Yang>(10, true);

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

        // Step 2: compute GRM = Z * Z^T
        Eigen::MatrixXd expected_grm = Z * Z.transpose();

        // Step 3: normalize by trace/n
        auto n = static_cast<double>(genotypes.rows());
        expected_grm /= expected_grm.trace() / n;

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
        Eigen::MatrixXd G = grm.compute<grm::Su>(10, true);

        // Su additive: mean centering
        Eigen::MatrixXd Z = genotypes;
        for (Eigen::Index j = 0; j < Z.cols(); ++j)
        {
            Z.col(j).array() -= Z.col(j).mean();
        }

        Eigen::MatrixXd expected_grm = Z * Z.transpose();
        auto n = static_cast<double>(genotypes.rows());
        expected_grm /= expected_grm.trace() / n;

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

        Eigen::MatrixXd G = grm.compute<grm::Yang>(5, true, callback);

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
    }

    SECTION("Null callback does not crash")
    {
        const Eigen::Index num_samples = 4;
        const Eigen::Index num_snps = 10;
        auto [bed_prefix, genotypes]
            = fixture.create_bed_files(num_samples, num_snps, 0.0);

        GRM grm(bed_prefix);

        REQUIRE_NOTHROW(grm.compute<grm::Yang>(5, true, nullptr));
    }
}
