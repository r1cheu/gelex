#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>

#include "../src/predictor/snp_effect_processor.h"

using namespace gelex;
using namespace Catch::Matchers;

class SnpEffectProcessorTestFixture
{
   public:
    SnpEffectProcessorTestFixture()
    {
        // Create a temporary .snp.eff file for testing
        test_snp_eff_file_ = std::filesystem::temp_directory_path()
                             / "test_snp_effects.snp.eff";

        std::ofstream file(test_snp_eff_file_);
        file << "Index\tID\tChrom\tPosition\tA1\tA2\tA1Frq\tAdd\tAddSE\tAddPVE"
                "\tPIP\tDom\tDomSE\tDomPVE\tPIP\n";
        file << "1\trs123\t1\t1000\tA\tT\t0.3\t0.5\t0.1\t0.01\t0.9\t0.2\t0."
                "05\t0.005\t0.8\n";
        file << "2\trs456\t1\t2000\tC\tG\t0.7\t-0.3\t0.08\t0.008\t0.95\t0.1\t0."
                "03\t0.003\t0.7\n";
        file << "3\trs789\t2\t3000\tG\tA\t0.5\t0.8\t0.12\t0.015\t0.98\t0.0\t0."
                "0\t0.0\t0.0\n";
        file << "4\trs101\t2\t4000\tT\tC\t0.2\t-0.6\t0.15\t0.012\t0."
                "92\tnan\tnan\tnan\tnan\n";
        file.close();
    }

    ~SnpEffectProcessorTestFixture()
    {
        // Clean up the temporary file
        if (std::filesystem::exists(test_snp_eff_file_))
        {
            std::filesystem::remove(test_snp_eff_file_);
        }
    }

    std::filesystem::path test_snp_eff_file_;
};

TEST_CASE(
    "SnpEffectProcessor calculates GEVI for single genotype",
    "[snp_effect_processor]")
{
    SnpEffect info{"rs123", "1", 1000, 'A', 'T', 0.3, 0.5, 0.2};

    SECTION("Genotype 0 - homozygous for A2")
    {
        double gevi = SnpEffectProcessor::calculate_gevi(0, info);

        // Manual calculation for verification:
        // p = 0.3, q = 0.7
        // add_encoded = 0, dom_encoded = 0
        // add_std = (0 - 2*0.3) / sqrt(2*0.3*0.7) = -0.6 / sqrt(0.42) ≈ -0.6 /
        // 0.648074 ≈ -0.92582 dom_std = (0 - 2*0.3*0.3) / (2*0.3*0.7) = (0 -
        // 0.18) / 0.42 = -0.42857 GEVI = (-0.92582 * 0.5) + (-0.42857 * 0.2) =
        // -0.46291 - 0.085714 ≈ -0.54862

        REQUIRE_THAT(gevi, WithinAbs(-0.54862, 1e-4));
    }

    SECTION("Genotype 1 - heterozygous")
    {
        double gevi = SnpEffectProcessor::calculate_gevi(1, info);

        // Manual calculation:
        // add_encoded = 1, dom_encoded = 2*0.3 = 0.6
        // add_std = (1 - 0.6) / sqrt(0.42) = 0.4 / 0.648074 ≈ 0.61721
        // dom_std = (0.6 - 0.18) / 0.42 = 0.42 / 0.42 = 1.0
        // GEVI = (0.61721 * 0.5) + (1.0 * 0.2) = 0.308605 + 0.2 = 0.508605

        REQUIRE_THAT(gevi, WithinAbs(0.50861, 1e-4));
    }

    SECTION("Genotype 2 - homozygous for A1")
    {
        double gevi = SnpEffectProcessor::calculate_gevi(2, info);

        // Manual calculation:
        // add_encoded = 2, dom_encoded = 4*0.3 - 2 = 1.2 - 2 = -0.8
        // add_std = (2 - 0.6) / sqrt(0.42) = 1.4 / 0.648074 ≈ 2.16025
        // dom_std = (-0.8 - 0.18) / 0.42 = -0.98 / 0.42 ≈ -2.33333
        // GEVI = (2.16025 * 0.5) + (-2.33333 * 0.2) = 1.080125 - 0.466666 ≈
        // 0.61346

        REQUIRE_THAT(gevi, WithinAbs(0.61346, 1e-4));
    }
}

TEST_CASE(
    "SnpEffectProcessor handles edge cases in GEVI calculation",
    "[snp_effect_processor]")
{
    SECTION("Zero frequency (monomorphic)")
    {
        SnpEffect info{"rs000", "1", 1000, 'A', 'T', 0.0, 0.5, 0.2};

        // Should handle zero frequency without division by zero
        double gevi = SnpEffectProcessor::calculate_gevi(1, info);
        // Implementation should handle this gracefully
        REQUIRE(!std::isnan(gevi));
        REQUIRE(!std::isinf(gevi));
    }

    SECTION("One frequency (monomorphic)")
    {
        SnpEffect info{"rs111", "1", 1000, 'A', 'T', 1.0, 0.5, 0.2};

        double gevi = SnpEffectProcessor::calculate_gevi(1, info);
        // Implementation should handle this gracefully
        REQUIRE(!std::isnan(gevi));
        REQUIRE(!std::isinf(gevi));
    }

    SECTION("Zero effects")
    {
        SnpEffect info{"rs000", "1", 1000, 'A', 'T', 0.5, 0.0, 0.0};

        double gevi = SnpEffectProcessor::calculate_gevi(1, info);
        REQUIRE_THAT(gevi, WithinAbs(0.0, 1e-10));
    }
}

TEST_CASE("SnpEffectProcessor calculates batch GEVI", "[snp_effect_processor]")
{
    SnpEffect info{"rs123", "1", 1000, 'A', 'T', 0.3, 0.5, 0.2};
    std::vector<int> genotypes{0, 1, 2, 1, 0};

    auto results = SnpEffectProcessor::calculate_gevi_batch(genotypes, info);

    REQUIRE(results.size() == genotypes.size());

    // Verify first genotype (homozygous A2)
    REQUIRE_THAT(results[0], WithinAbs(-0.54862, 1e-4));

    // Verify second genotype (heterozygous)
    REQUIRE_THAT(results[1], WithinAbs(0.50861, 1e-4));

    // Verify third genotype (homozygous A1)
    REQUIRE_THAT(results[2], WithinAbs(0.61346, 1e-4));
}

TEST_CASE_METHOD(
    SnpEffectProcessorTestFixture,
    "SnpEffectProcessor loads .snp.eff file",
    "[snp_effect_processor]")
{
    auto snp_processor = SnpEffectProcessor::create(test_snp_eff_file_);
    REQUIRE(snp_processor.has_value());

    const auto& infos = snp_processor->snp_effects();
    REQUIRE(infos.size() == 4);

    // Verify first SNP
    REQUIRE(infos[0].meta.id == "rs123");
    REQUIRE(infos[0].meta.a1 == 'A');
    REQUIRE(infos[0].meta.a2 == 'T');
    REQUIRE_THAT(infos[0].p_freq, WithinAbs(0.3, 1e-6));
    REQUIRE_THAT(infos[0].add_effect, WithinAbs(0.5, 1e-6));
    REQUIRE_THAT(infos[0].dom_effect, WithinAbs(0.2, 1e-6));

    // Verify second SNP
    REQUIRE(infos[1].meta.id == "rs456");
    REQUIRE(infos[1].meta.a1 == 'C');
    REQUIRE(infos[1].meta.a2 == 'G');
    REQUIRE_THAT(infos[1].p_freq, WithinAbs(0.7, 1e-6));
    REQUIRE_THAT(infos[1].add_effect, WithinAbs(-0.3, 1e-6));
    REQUIRE_THAT(infos[1].dom_effect, WithinAbs(0.1, 1e-6));

    // Verify third SNP (dominant effect is zero)
    REQUIRE(infos[2].meta.id == "rs789");
    REQUIRE_THAT(infos[2].dom_effect, WithinAbs(0.0, 1e-6));

    // Verify fourth SNP (dominant effect not available, should be zero)
    REQUIRE(infos[3].meta.id == "rs101");
    REQUIRE(std::isnan(infos[3].dom_effect));
}

TEST_CASE("SnpEffectProcessor handles empty inputs", "[snp_effect_processor]")
{
    SECTION("Empty genotype vector")
    {
        std::vector<int> empty_genotypes;
        SnpEffect info{"rs123", "1", 1000, 'A', 'T', 0.5, 0.5, 0.2};

        auto results
            = SnpEffectProcessor::calculate_gevi_batch(empty_genotypes, info);
        REQUIRE(results.empty());
    }

    SECTION("Empty genotype matrix")
    {
        std::vector<std::vector<int>> empty_genotypes;
        std::vector<SnpEffect> empty_snp_infos;

        auto results = SnpEffectProcessor::calculate_total_genetic_value(
            empty_genotypes, empty_snp_infos);
        REQUIRE(results.empty());
    }

    SECTION("Mismatched dimensions")
    {
        std::vector<std::vector<int>> genotypes{
            {0, 1}, {1, 0}};                   // 2 SNPs, 2 individuals
        std::vector<SnpEffect> snp_infos{{}};  // Only 1 SNP

        auto results = SnpEffectProcessor::calculate_total_genetic_value(
            genotypes, snp_infos);
        REQUIRE(results.empty());
    }
}

TEST_CASE(
    "SnpEffectProcessor returns error for non-existent file",
    "[snp_effect_processor]")
{
    auto snp_infos = SnpEffectProcessor::create("non_existent_file.snp.eff");
    REQUIRE(!snp_infos.has_value());
}

TEST_CASE(
    "SnpEffectProcessor handles missing required columns",
    "[snp_effect_processor]")
{
    // Create a file with missing required columns
    std::filesystem::path missing_columns_file
        = std::filesystem::temp_directory_path() / "missing_columns.snp.eff";

    {
        std::ofstream file(missing_columns_file);
        file << "Index\tChrom\tPosition\tA1\tA2\tA1Frq\n";  // Missing ID and
                                                            // Add columns
        file << "1\t1\t1000\tA\tT\t0.3\n";
        file.close();
    }

    auto snp_infos = SnpEffectProcessor::create(missing_columns_file);
    REQUIRE(!snp_infos.has_value());

    // Clean up
    std::filesystem::remove(missing_columns_file);
}
