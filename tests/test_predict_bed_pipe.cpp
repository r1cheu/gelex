#include <fstream>
#include <sstream>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../src/predictor/predict_bed_pipe.h"
#include "../src/predictor/snp_matcher.h"
#include "bed_fixture.h"
#include "file_fixture.h"
#include "gelex/data/sample_manager.h"
#include "gelex/exception.h"

namespace fs = std::filesystem;

using namespace gelex;  // For PredictBedPipe, MatchType, etc.
using Catch::Matchers::EndsWith;
using gelex::test::are_matrices_equal;
using gelex::test::BedFixture;
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

TEST_CASE("PredictBedPipe - Constructor", "[predictor][predict_bed_pipe]")
{
    SECTION("Happy path - successful construction with valid files")
    {
        // Create BED files using BedFixture
        BedFixture bed_fixture;
        const Eigen::Index num_samples = 3;
        const Eigen::Index num_snps = 3;
        auto [bed_prefix, _]
            = bed_fixture.create_bed_files(num_samples, num_snps, 0.0);

        // Create SNP effect file in the same directory as BED files
        auto& file_fixture = bed_fixture.get_file_fixture();
        auto snp_effects = create_snp_effects(
            file_fixture,
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs001\tA\tC\t0.25\t0.123\t0.045",
             "rs002\tT\tG\t0.75\t-0.456\t0.089",
             "rs003\tC\tA\t0.50\t0.789\t-0.012"});

        // Update BIM file to match SNP effect file
        auto bim_path = bed_prefix;
        bim_path.replace_extension(".bim");
        {
            std::ofstream bim_file(bim_path);
            bim_file << "1\trs001\t0\t1000\tA\tC\n";
            bim_file << "1\trs002\t0\t2000\tT\tG\n";
            bim_file << "1\trs003\t0\t3000\tC\tA\n";
        }

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto sample_manager = std::make_shared<SampleManager>(fam_path, false);
        sample_manager->finalize();

        REQUIRE_NOTHROW(
            [&]()
            {
                PredictBedPipe pipe(bed_prefix, snp_effects, sample_manager);
            }());
    }

    SECTION("Exception path - sample manager is nullptr")
    {
        BedFixture bed_fixture;
        const Eigen::Index num_samples = 3;
        const Eigen::Index num_snps = 3;
        auto [bed_prefix, _]
            = bed_fixture.create_bed_files(num_samples, num_snps, 0.0);

        // Create SNP effect file in the same directory as BED files
        auto& file_fixture = bed_fixture.get_file_fixture();
        auto snp_effects = create_snp_effects(
            file_fixture,
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs001\tA\tC\t0.25\t0.123\t0.045"});

        REQUIRE_THROWS_MATCHES(
            PredictBedPipe(bed_prefix, snp_effects, nullptr),
            ArgumentValidationException,
            Catch::Matchers::MessageMatches(
                EndsWith("SampleManager cannot be null")));
    }
}

TEST_CASE(
    "PredictBedPipe - load() method with MatchPlan filtering",
    "[predictor][predict_bed_pipe]")
{
    // Core test for MatchPlan SNP filtering correctness

    SECTION("Scenario A: Perfect match (all keep)")
    {
        // Create BED file with 3 samples and 3 SNPs
        BedFixture bed_fixture;
        const Eigen::Index num_samples = 3;
        const Eigen::Index num_snps = 3;

        auto [bed_prefix, genotypes]
            = bed_fixture.create_bed_files(num_samples, num_snps, 0.0);

        // Create SNP effect file matching all SNPs with same alleles
        // Use the same FileFixture as BedFixture to keep files in same
        // directory
        auto& file_fixture = bed_fixture.get_file_fixture();
        auto snp_effects = create_snp_effects(
            file_fixture,
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs001\tA\tC\t0.25\t0.123\t0.045",
             "rs002\tT\tG\t0.75\t-0.456\t0.089",
             "rs003\tC\tA\t0.50\t0.789\t-0.012"});

        // Update BIM file to match alleles (BedFixture creates random SNP IDs,
        // we need to update them)
        auto bim_path = bed_prefix;
        bim_path.replace_extension(".bim");
        {
            std::ofstream bim_file(bim_path);
            bim_file << "1\trs001\t0\t1000\tA\tC\n";
            bim_file << "1\trs002\t0\t2000\tT\tG\n";
            bim_file << "1\trs003\t0\t3000\tC\tA\n";
        }

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        REQUIRE(fs::exists(fam_path));
        auto sample_manager = std::make_shared<SampleManager>(fam_path, false);
        sample_manager->finalize();

        PredictBedPipe pipe(bed_prefix, snp_effects, sample_manager);

        // Load the filtered matrix
        Eigen::MatrixXd filtered = pipe.load();

        // Verify dimensions
        REQUIRE(filtered.rows() == num_samples);
        REQUIRE(
            filtered.cols()
            == num_snps);  // Should have same number of columns as effect SNPs
        are_matrices_equal(filtered, genotypes, 1e-8);
    }

    SECTION("Scenario B: Reverse match")
    {
        // Test reverse matching: col = 2.0 - col.array()

        BedFixture bed_fixture;
        const Eigen::Index num_samples = 2;
        const Eigen::Index num_snps = 1;

        auto [bed_prefix, genotypes]
            = bed_fixture.create_bed_files(num_samples, num_snps, 0.0);

        // Load original genotypes using BedPipe
        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        REQUIRE(fs::exists(fam_path));
        auto sample_manager = std::make_shared<SampleManager>(fam_path, false);
        sample_manager->finalize();

        BedPipe bed_pipe(bed_prefix, sample_manager);
        Eigen::MatrixXd original = bed_pipe.load();

        // Create SNP effect file in the same directory as BED files
        auto& file_fixture = bed_fixture.get_file_fixture();
        auto snp_effects = create_snp_effects(
            file_fixture,
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs001\tA\tC\t0.25\t0.123\t0.045"});

        // Update BIM file with swapped alleles (to trigger reverse match)
        auto bim_path = bed_prefix;
        bim_path.replace_extension(".bim");
        {
            std::ofstream bim_file(bim_path);
            bim_file << "1\trs001\t0\t1000\tC\tA\n";  // Swapped!
        }

        PredictBedPipe pipe(bed_prefix, snp_effects, sample_manager);
        Eigen::MatrixXd filtered = pipe.load();

        // Calculate expected reverse: 2.0 - original
        Eigen::MatrixXd expected = 2.0 - original.array();

        // Verify dimensions
        REQUIRE(filtered.rows() == num_samples);
        REQUIRE(filtered.cols() == 1);
        REQUIRE(are_matrices_equal(filtered, expected, 1e-8));
    }

    SECTION("Scenario C: Mixed match types (keep, reverse, skip)")
    {
        // Test a mixture of match types

        BedFixture bed_fixture;
        const Eigen::Index num_samples = 2;
        const Eigen::Index num_snps = 3;

        auto [bed_prefix, _]
            = bed_fixture.create_bed_files(num_samples, num_snps, 0.0);

        // Load original genotypes
        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        REQUIRE(fs::exists(fam_path));
        auto sample_manager = std::make_shared<SampleManager>(fam_path, false);
        sample_manager->finalize();

        BedPipe bed_pipe(bed_prefix, sample_manager);
        Eigen::MatrixXd original = bed_pipe.load();

        // Create SNP effect file with 2 SNPs (one will be missing)
        auto& file_fixture = bed_fixture.get_file_fixture();
        auto snp_effects = create_snp_effects(
            file_fixture,
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs001\tA\tC\t0.25\t0.123\t0.045",     // Will match with keep
             "rs002\tT\tG\t0.75\t-0.456\t0.089"});  // Will match with reverse

        // Update BIM file with mixed scenarios:
        // SNP1: rs001 A C (keep)
        // SNP2: rs002 G T (reverse - alleles swapped)
        // SNP3: rs003 A G (skip - SNP not in effects)
        auto bim_path = bed_prefix;
        bim_path.replace_extension(".bim");
        {
            std::ofstream bim_file(bim_path);
            bim_file << "1\trs001\t0\t1000\tA\tC\n";  // keep
            bim_file
                << "1\trs002\t0\t2000\tG\tT\n";  // reverse (T/G swapped to G/T)
            bim_file << "1\trs003\t0\t3000\tA\tG\n";  // skip
        }

        PredictBedPipe pipe(bed_prefix, snp_effects, sample_manager);
        Eigen::MatrixXd filtered = pipe.load();

        // Verify dimensions: should have 2 columns (2 effect SNPs)
        REQUIRE(filtered.rows() == num_samples);
        REQUIRE(filtered.cols() == 2);

        // Check that column 0 (rs001) matches original column 0 (keep)
        for (Eigen::Index i = 0; i < num_samples; ++i)
        {
            double orig = original(i, 0);
            double filt = filtered(i, 0);
            REQUIRE(filt == orig);
        }

        // Check that column 1 (rs002) is reverse of original column 1
        for (Eigen::Index i = 0; i < num_samples; ++i)
        {
            double orig = original(i, 1);
            double filt = filtered(i, 1);
            REQUIRE(filt == 2.0 - orig);
        }
    }

    SECTION("Scenario D: No matching SNPs (all skip)")
    {
        // Test when no SNPs match the effect file

        BedFixture bed_fixture;
        const Eigen::Index num_samples = 3;
        const Eigen::Index num_snps = 2;

        auto [bed_prefix, _]
            = bed_fixture.create_bed_files(num_samples, num_snps, 0.0);

        // Create SNP effect file with different SNPs
        auto& file_fixture = bed_fixture.get_file_fixture();
        auto snp_effects = create_snp_effects(
            file_fixture,
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs999\tA\tC\t0.25\t0.123\t0.045",  // Different SNP ID
             "rs998\tT\tG\t0.75\t-0.456\t0.089"});

        // Update BIM file with different SNPs
        auto bim_path = bed_prefix;
        bim_path.replace_extension(".bim");
        {
            std::ofstream bim_file(bim_path);
            bim_file << "1\trs001\t0\t1000\tA\tC\n";  // Not in effects
            bim_file << "1\trs002\t0\t2000\tT\tG\n";  // Not in effects
        }

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        REQUIRE(fs::exists(fam_path));
        auto sample_manager = std::make_shared<SampleManager>(fam_path, false);
        sample_manager->finalize();

        PredictBedPipe pipe(bed_prefix, snp_effects, sample_manager);
        Eigen::MatrixXd filtered = pipe.load();

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

TEST_CASE("PredictBedPipe - Edge cases", "[predictor][predict_bed_pipe]")
{
    SECTION("Sparse SampleManager interaction")
    {
        // Test with SampleManager that filters samples

        BedFixture bed_fixture;
        const Eigen::Index num_raw_samples = 5;
        const Eigen::Index num_snps = 4;

        auto [bed_prefix, genotypes]
            = bed_fixture.create_bed_files(num_raw_samples, num_snps, 0.0);

        // Create SampleManager that only includes first 2 samples
        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        REQUIRE(fs::exists(fam_path));
        auto sample_manager = std::make_shared<SampleManager>(fam_path, true);

        // Read sample IDs from FAM file
        std::ifstream fam_file(fam_path);
        std::vector<std::string> raw_ids;
        std::string line;
        while (std::getline(fam_file, line))
        {
            std::istringstream iss(line);
            std::string iid;
            iss >> iid >> iid;
            raw_ids.push_back(iid);
        }

        // Intersect with first 2 samples
        std::vector<std::string_view> intersect_ids;
        for (size_t i = 0; i < 2; ++i)
        {
            intersect_ids.push_back(raw_ids[i]);
        }
        sample_manager->intersect(intersect_ids);
        sample_manager->finalize();

        // Create SNP effect file matching all SNPs
        auto& file_fixture = bed_fixture.get_file_fixture();
        auto snp_effects = create_snp_effects(
            file_fixture,
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs001\tA\tC\t0.25\t0.123\t0.045",
             "rs002\tT\tG\t0.75\t-0.456\t0.089",
             "rs003\tC\tA\t0.50\t0.789\t-0.012"});

        // Update BIM file with matching alleles
        auto bim_path = bed_prefix;
        bim_path.replace_extension(".bim");
        {
            std::ofstream bim_file(bim_path);
            bim_file << "1\trs001\t0\t1000\ta\tc\n";  // lower case alleles
            bim_file
                << "1\trs002\t0\t2000\tg\tT\n";  // one lower case with reverse
            bim_file << "1\trs003\t0\t3000\tC\tG\n";  // skipped SNP
            bim_file << "1\trs004\t0\t4000\tC\tG\n";  // skipped SNP
        }

        PredictBedPipe pipe(bed_prefix, snp_effects, sample_manager);
        Eigen::MatrixXd filtered = pipe.load();
        genotypes
            = genotypes.block(0, 0, 2, 3).eval();  // Keep only first 2 samples
        genotypes.col(1)
            = 2 - genotypes.col(1).array();  // Reverse match for SNP2
        genotypes.col(2).setZero();

        // Verify dimensions: should have 2 rows (filtered samples) and 3
        // columns
        REQUIRE(filtered.rows() == 2);
        REQUIRE(filtered.cols() == 3);
        REQUIRE(are_matrices_equal(filtered, genotypes, 1e-8));
    }
}
