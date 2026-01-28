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

#include <fstream>
#include <sstream>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../src/data/loader/snp_effect_loader.h"
#include "../src/predict/genotype_aligner.h"
#include "../src/predict/snp_matcher.h"

#include "bed_fixture.h"
#include "file_fixture.h"

#include "gelex/exception.h"

namespace fs = std::filesystem;

using namespace gelex;  // For GenotypeAligner, MatchType, etc.
using Catch::Matchers::EndsWith;
using gelex::test::are_matrices_equal;

using gelex::test::FileFixture;

namespace
{

// Helper function to create SNP effect test files (reused from
// test_snp_matcher.cpp)
std::string create_snp_effect_content(
    const std::string& header,
    const std::vector<std::string>& rows)
{
    std::string content = header + "\n";
    for (const auto& row : rows)
    {
        content += row + "\n";
    }
    return content;
}

// Helper function to create SnpEffects from content
static SnpEffects create_snp_effects(
    FileFixture& files,
    const std::string& header,
    const std::vector<std::string>& rows)
{
    std::string content = create_snp_effect_content(header, rows);
    auto file_path = files.create_text_file(content, ".snp.eff");
    gelex::detail::SnpEffectLoader loader(file_path);
    return std::move(loader).take_effects();
}

}  // namespace

TEST_CASE("GenotypeAligner - Constructor", "[predict][genotype_aligner]")
{
    SECTION("Happy path - successful construction with valid files")
    {
        FileFixture file_fixture;

        // Create .bim file with SNP information
        std::string bim_content
            = "1\trs001\t0\t1000\tA\tC\n"
              "1\trs002\t0\t2000\tT\tG\n"
              "1\trs003\t0\t3000\tC\tA\n";
        auto bim_path = file_fixture.create_text_file(bim_content, ".bim");

        // Create matching SNP effects
        auto snp_effects = create_snp_effects(
            file_fixture,
            "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
            {"1\t1000\trs001\tA\tC\t0.25\t0.123\t0.045",
             "1\t2000\trs002\tT\tG\t0.75\t-0.456\t0.089",
             "1\t3000\trs003\tC\tA\t0.50\t0.789\t-0.012"});

        // Construct GenotypeAligner - should not throw
        REQUIRE_NOTHROW(
            [&]()
            {
                // Use bim path without .bim extension as bed prefix
                auto bed_path = bim_path;
                bed_path.replace_extension("bed");
                GenotypeAligner aligner(bed_path, snp_effects);
            }());
    }
}

TEST_CASE(
    "GenotypeAligner - load() method with MatchPlan filtering",
    "[predict][genotype_aligner]")
{
    // Core test for MatchPlan SNP filtering correctness

    SECTION("Scenario A: Perfect match (all keep)")
    {
        FileFixture file_fixture;
        const Eigen::Index num_samples = 3;
        const Eigen::Index num_snps = 3;

        // Create .bim file with 3 SNPs
        std::string bim_content
            = "1\trs001\t0\t1000\tA\tC\n"
              "1\trs002\t0\t2000\tT\tG\n"
              "1\trs003\t0\t3000\tC\tA\n";
        auto bim_path = file_fixture.create_text_file(bim_content, ".bim");
        auto bed_prefix = bim_path;
        bed_prefix.replace_extension("");  // Remove .bim extension

        // Create matching SNP effects
        auto snp_effects = create_snp_effects(
            file_fixture,
            "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
            {"1\t1000\trs001\tA\tC\t0.25\t0.123\t0.045",
             "1\t2000\trs002\tT\tG\t0.75\t-0.456\t0.089",
             "1\t3000\trs003\tC\tA\t0.50\t0.789\t-0.012"});

        // Create genotype matrix with fixed values
        Eigen::MatrixXd genotypes(num_samples, num_snps);
        genotypes << 0.0, 1.0, 2.0, 1.0, 2.0, 0.0, 2.0, 0.0, 1.0;

        // Create GenotypeAligner and load aligned genotypes
        GenotypeAligner aligner(bed_prefix, snp_effects);
        Eigen::MatrixXd original_genotypes
            = genotypes;  // Keep copy for comparison
        Eigen::MatrixXd filtered = aligner.align(std::move(genotypes));

        // Verify dimensions
        REQUIRE(filtered.rows() == num_samples);
        REQUIRE(
            filtered.cols()
            == num_snps);  // Should have same number of columns as effect SNPs

        // Verify all SNPs are kept (perfect match)
        REQUIRE(are_matrices_equal(filtered, original_genotypes, 1e-8));
    }

    SECTION("Scenario B: Reverse match")
    {
        // Test reverse matching: col = 2.0 - col.array()
        FileFixture file_fixture;
        const Eigen::Index num_samples = 2;
        const Eigen::Index num_snps = 1;

        // Create .bim file with swapped alleles (C A instead of A C)
        std::string bim_content = "1\trs001\t0\t1000\tC\tA\n";  // Swapped!
        auto bim_path = file_fixture.create_text_file(bim_content, ".bim");
        auto bed_path = bim_path;
        bed_path.replace_extension("bed");  // Remove .bim extension

        // Create SNP effect file with original alleles (A C)
        auto snp_effects = create_snp_effects(
            file_fixture,
            "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
            {"1\t1000\trs001\tA\tC\t0.25\t0.123\t0.045"});

        // Create genotype matrix with fixed values
        Eigen::MatrixXd genotypes(num_samples, num_snps);
        genotypes << 0.0,  // Sample 1
            2.0;           // Sample 2

        // Calculate expected reverse: 2.0 - original
        Eigen::MatrixXd expected = 2.0 - genotypes.array();

        // Create GenotypeAligner and load aligned genotypes
        GenotypeAligner aligner(bed_path, snp_effects);
        Eigen::MatrixXd filtered = aligner.align(std::move(genotypes));

        // Verify dimensions
        REQUIRE(filtered.rows() == num_samples);
        REQUIRE(filtered.cols() == 1);
        REQUIRE(are_matrices_equal(filtered, expected, 1e-8));
    }

    SECTION("Scenario C: Mixed match types (keep, reverse, skip)")
    {
        // Test a mixture of match types
        FileFixture file_fixture;
        const Eigen::Index num_samples = 2;
        const Eigen::Index num_snps = 3;

        // Create .bim file with mixed scenarios:
        // SNP1: rs001 A C (keep)
        // SNP2: rs002 G T (reverse - alleles swapped)
        // SNP3: rs003 A G (skip - SNP not in effects)
        std::string bim_content
            = "1\trs001\t0\t1000\tA\tC\n"   // keep
              "1\trs002\t0\t2000\tG\tT\n"   // reverse (T/G swapped to G/T)
              "1\trs003\t0\t3000\tA\tG\n";  // skip
        auto bim_path = file_fixture.create_text_file(bim_content, ".bim");
        auto bed_prefix = bim_path;
        bed_prefix.replace_extension("");  // Remove .bim extension

        // Create SNP effect file with 2 SNPs (one will be missing)
        auto snp_effects = create_snp_effects(
            file_fixture,
            "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
            {"1\t1000\trs001\tA\tC\t0.25\t0.123\t0.045",     // Will match with
                                                             // keep
             "1\t2000\trs002\tT\tG\t0.75\t-0.456\t0.089"});  // Will match with
                                                             // reverse

        // Create genotype matrix with fixed values
        Eigen::MatrixXd genotypes(num_samples, num_snps);
        genotypes << 0.0, 1.0, 2.0,  // Sample 1
            1.0, 2.0, 0.0;           // Sample 2

        // Create GenotypeAligner and load aligned genotypes
        GenotypeAligner aligner(bed_prefix, snp_effects);
        Eigen::MatrixXd original_genotypes
            = genotypes;  // Keep copy for comparison
        Eigen::MatrixXd filtered = aligner.align(std::move(genotypes));

        // Verify dimensions: should have 2 columns (2 effect SNPs)
        REQUIRE(filtered.rows() == num_samples);
        REQUIRE(filtered.cols() == 2);

        // Check that column 0 (rs001) matches original column 0 (keep)
        for (Eigen::Index i = 0; i < num_samples; ++i)
        {
            double orig = original_genotypes(i, 0);
            double filt = filtered(i, 0);
            REQUIRE(filt == orig);
        }

        // Check that column 1 (rs002) is reverse of original column 1
        for (Eigen::Index i = 0; i < num_samples; ++i)
        {
            double orig = original_genotypes(i, 1);
            double filt = filtered(i, 1);
            REQUIRE(filt == 2.0 - orig);
        }
    }

    SECTION("Scenario D: No matching SNPs (all skip)")
    {
        // Test when no SNPs match the effect file
        FileFixture file_fixture;
        const Eigen::Index num_samples = 3;
        const Eigen::Index num_snps = 2;

        // Create .bim file with SNPs not in effects
        std::string bim_content
            = "1\trs001\t0\t1000\tA\tC\n"   // Not in effects
              "1\trs002\t0\t2000\tT\tG\n";  // Not in effects
        auto bim_path = file_fixture.create_text_file(bim_content, ".bim");
        auto bed_prefix = bim_path;
        bed_prefix.replace_extension("");  // Remove .bim extension

        // Create SNP effect file with different SNPs
        auto snp_effects = create_snp_effects(
            file_fixture,
            "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
            {"1\t1000\trs999\tA\tC\t0.25\t0.123\t0.045",  // Different SNP ID
             "1\t2000\trs998\tT\tG\t0.75\t-0.456\t0.089"});

        // Create genotype matrix with fixed values
        Eigen::MatrixXd genotypes(num_samples, num_snps);
        genotypes << 0.0, 1.0, 1.0, 2.0, 2.0, 0.0;

        // Create GenotypeAligner and load aligned genotypes
        GenotypeAligner aligner(bed_prefix, snp_effects);
        Eigen::MatrixXd filtered = aligner.align(std::move(genotypes));

        // Verify dimensions: should have 2 columns (2 effect SNPs)
        REQUIRE(filtered.rows() == num_samples);
        REQUIRE(filtered.cols() == 2);

        // All values should be zero (since all SNPs are skipped and matrix was
        // setZero)
        for (Eigen::Index j = 0; j < filtered.cols(); ++j)
        {
            for (Eigen::Index i = 0; i < filtered.rows(); ++i)
            {
                REQUIRE(filtered(i, j) == 0.0);
            }
        }
    }
}
