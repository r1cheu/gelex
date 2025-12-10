#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "file_fixture.h"

#include "../src/predictor/snp_matcher.h"

namespace fs = std::filesystem;

using namespace gelex;  // For SnpMatcher, MatchType, etc.
using Catch::Matchers::EndsWith;
using gelex::test::FileFixture;

// Helper function to create SNP effect test files
static std::string create_snp_effect_content(
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

// Helper function to create BIM test files
static std::string create_bim_content(const std::vector<std::string>& rows)
{
    std::string content;
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

TEST_CASE("SnpMatcher - Constructor", "[predictor][snp_matcher]")
{
    FileFixture files;

    SECTION("Happy path - successful construction with valid SNP effect file")
    {
        auto effects = create_snp_effects(
            files,
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs001\tA\tC\t0.25\t0.123\t0.045",
             "rs002\tT\tG\t0.75\t-0.456\t0.089",
             "rs003\tC\tA\t0.50\t0.789\t-0.012"});

        REQUIRE_NOTHROW([&]() { SnpMatcher matcher(effects); }());
    }
}

TEST_CASE("SnpMatcher - match() method", "[predictor][snp_matcher]")
{
    FileFixture files;

    SECTION("Happy path - perfect match (all alleles identical)")
    {
        // Create SNP effects
        auto effects = create_snp_effects(
            files,
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs001\tA\tC\t0.25\t0.123\t0.045",
             "rs002\tT\tG\t0.75\t-0.456\t0.089",
             "rs003\tC\tA\t0.50\t0.789\t-0.012"});

        // Create BIM file with same SNPs and alleles
        std::string bim_content = create_bim_content(
            {"1\trs001\t0\t1000\tA\tC",
             "1\trs002\t0\t2000\tT\tG",
             "1\trs003\t0\t3000\tC\tA"});

        auto bim_path = files.create_named_text_file("test.bim", bim_content);
        auto bed_path = bim_path;
        bed_path.replace_extension(".bed");

        SnpMatcher matcher(effects);
        auto match_plan = matcher.match(bed_path);

        REQUIRE(match_plan.size() == 3);

        // Check all matches are "keep" with correct target columns
        REQUIRE(match_plan[0].type == MatchType::keep);
        REQUIRE(
            match_plan[0].target_col == 0);  // rs001 is first in SNP effects

        REQUIRE(match_plan[1].type == MatchType::keep);
        REQUIRE(match_plan[1].target_col == 1);  // rs002 is second

        REQUIRE(match_plan[2].type == MatchType::keep);
        REQUIRE(match_plan[2].target_col == 2);  // rs003 is third
    }

    SECTION("Happy path - reverse match (alleles swapped)")
    {
        // Create SNP effects
        auto effects = create_snp_effects(
            files,
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs001\tA\tC\t0.25\t0.123\t0.045",
             "rs002\tT\tG\t0.75\t-0.456\t0.089"});

        // Create BIM file with swapped alleles
        std::string bim_content = create_bim_content(
            {"1\trs001\t0\t1000\tC\tA",    // A/C swapped to C/A
             "1\trs002\t0\t2000\tG\tT"});  // T/G swapped to G/T

        auto bim_path = files.create_named_text_file("test.bim", bim_content);
        auto bed_path = bim_path;
        bed_path.replace_extension(".bed");

        SnpMatcher matcher(effects);
        auto match_plan = matcher.match(bed_path);

        REQUIRE(match_plan.size() == 2);

        // Check all matches are "reverse"
        REQUIRE(match_plan[0].type == MatchType::reverse);
        REQUIRE(match_plan[0].target_col == 0);

        REQUIRE(match_plan[1].type == MatchType::reverse);
        REQUIRE(match_plan[1].target_col == 1);
    }

    SECTION("Happy path - partial match (some match, some skip)")
    {
        // Create SNP effects
        auto effects = create_snp_effects(
            files,
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs001\tA\tC\t0.25\t0.123\t0.045",
             "rs002\tT\tG\t0.75\t-0.456\t0.089",
             "rs003\tC\tA\t0.50\t0.789\t-0.012"});

        // Create BIM file with mixed matches
        std::string bim_content = create_bim_content(
            {"1\trs001\t0\t1000\tA\tC",    // keep
             "1\trs002\t0\t2000\tG\tT",    // reverse
             "1\trs003\t0\t3000\tA\tG",    // skip (different alleles)
             "1\trs004\t0\t4000\tT\tC"});  // skip (SNP not in effects)

        auto bim_path = files.create_named_text_file("test.bim", bim_content);
        auto bed_path = bim_path;
        bed_path.replace_extension(".bed");

        SnpMatcher matcher(effects);
        auto match_plan = matcher.match(bed_path);

        REQUIRE(match_plan.size() == 4);

        // Check match types
        REQUIRE(match_plan[0].type == MatchType::keep);
        REQUIRE(match_plan[0].target_col == 0);

        REQUIRE(match_plan[1].type == MatchType::reverse);
        REQUIRE(match_plan[1].target_col == 1);

        REQUIRE(match_plan[2].type == MatchType::skip);
        REQUIRE(match_plan[2].target_col == -1);

        REQUIRE(match_plan[3].type == MatchType::skip);
        REQUIRE(match_plan[3].target_col == -1);
    }

    SECTION("Happy path - case insensitive allele matching")
    {
        // Create SNP effects with uppercase alleles
        auto effects = create_snp_effects(
            files,
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs001\tA\tC\t0.25\t0.123\t0.045",
             "rs002\tT\tG\t0.75\t-0.456\t0.089"});

        // Create BIM file with lowercase alleles
        std::string bim_content = create_bim_content(
            {"1\trs001\t0\t1000\ta\tc",    // lowercase
             "1\trs002\t0\t2000\tt\tg"});  // lowercase

        auto bim_path = files.create_named_text_file("test.bim", bim_content);
        auto bed_path = bim_path;
        bed_path.replace_extension(".bed");

        SnpMatcher matcher(effects);
        auto match_plan = matcher.match(bed_path);

        REQUIRE(match_plan.size() == 2);

        // Check matches are "keep" despite case difference
        REQUIRE(match_plan[0].type == MatchType::keep);
        REQUIRE(match_plan[0].target_col == 0);

        REQUIRE(match_plan[1].type == MatchType::keep);
        REQUIRE(match_plan[1].target_col == 1);
    }

    SECTION("Happy path - no matching SNPs")
    {
        // Create SNP effects
        auto effects = create_snp_effects(
            files,
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs001\tA\tC\t0.25\t0.123\t0.045",
             "rs002\tT\tG\t0.75\t-0.456\t0.089"});

        // Create BIM file with completely different SNPs
        std::string bim_content = create_bim_content(
            {"1\trs101\t0\t1000\tA\tC",
             "1\trs102\t0\t2000\tT\tG",
             "1\trs103\t0\t3000\tC\tA"});

        auto bim_path = files.create_named_text_file("test.bim", bim_content);
        auto bed_path = bim_path;
        bed_path.replace_extension(".bed");

        SnpMatcher matcher(effects);
        auto match_plan = matcher.match(bed_path);

        REQUIRE(match_plan.size() == 3);

        // All should be skip
        for (const auto& info : match_plan.plan)
        {
            REQUIRE(info.type == MatchType::skip);
            REQUIRE(info.target_col == -1);
        }
    }
}

TEST_CASE(
    "SnpMatcher - determine_match_type() logic",
    "[predictor][snp_matcher]")
{
    // Note: Since determine_match_type() is private, we test it indirectly
    // through the public match() method

    FileFixture files;

    SECTION("Test allele combinations - keep cases")
    {
        // Create SNP effects
        auto effects = create_snp_effects(
            files,
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs001\tA\tC\t0.25\t0.123\t0.045",
             "rs002\ta\tc\t0.25\t0.123\t0.045"});  // lowercase

        // Test various keep scenarios
        std::string bim_content = create_bim_content(
            {"1\trs001\t0\t1000\tA\tC",    // uppercase match
             "1\trs002\t0\t2000\ta\tc"});  // lowercase match

        auto bim_path = files.create_named_text_file("test.bim", bim_content);
        auto bed_path = bim_path;
        bed_path.replace_extension(".bed");

        SnpMatcher matcher(effects);
        auto match_plan = matcher.match(bed_path);

        REQUIRE(match_plan.size() == 2);
        REQUIRE(match_plan[0].type == MatchType::keep);
        REQUIRE(match_plan[1].type == MatchType::keep);
    }

    SECTION("Test allele combinations - reverse cases")
    {
        // Create SNP effects
        auto effects = create_snp_effects(
            files,
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs001\tA\tC\t0.25\t0.123\t0.045",
             "rs002\ta\tc\t0.25\t0.123\t0.045"});  // lowercase

        // Test various reverse scenarios
        std::string bim_content = create_bim_content(
            {"1\trs001\t0\t1000\tC\tA",    // uppercase reverse
             "1\trs002\t0\t2000\tc\ta"});  // lowercase reverse

        auto bim_path = files.create_named_text_file("test.bim", bim_content);
        auto bed_path = bim_path;
        bed_path.replace_extension(".bed");

        SnpMatcher matcher(effects);
        auto match_plan = matcher.match(bed_path);

        REQUIRE(match_plan.size() == 2);
        REQUIRE(match_plan[0].type == MatchType::reverse);
        REQUIRE(match_plan[1].type == MatchType::reverse);
    }

    SECTION("Test allele combinations - skip cases")
    {
        // Create SNP effects
        auto effects = create_snp_effects(
            files,
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs001\tA\tC\t0.25\t0.123\t0.045"});

        // Test various skip scenarios
        std::string bim_content = create_bim_content(
            {"1\trs001\t0\t1000\tA\tG",    // A2 different
             "1\trs001\t0\t2000\tT\tC",    // A1 different
             "1\trs001\t0\t3000\tT\tG"});  // both different

        auto bim_path = files.create_named_text_file("test.bim", bim_content);
        auto bed_path = bim_path;
        bed_path.replace_extension(".bed");

        SnpMatcher matcher(effects);
        auto match_plan = matcher.match(bed_path);

        REQUIRE(match_plan.size() == 3);
        for (const auto& info : match_plan.plan)
        {
            REQUIRE(info.type == MatchType::skip);
            REQUIRE(info.target_col == -1);
        }
    }

    SECTION("Test allele combinations - case mixing")
    {
        // Create SNP effects with mixed case
        auto effects = create_snp_effects(
            files,
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs001\tA\tc\t0.25\t0.123\t0.045"});  // A uppercase, c lowercase

        // Test case mixing scenarios
        std::string bim_content = create_bim_content(
            {"1\trs001\t0\t1000\ta\tC",    // a lowercase, C uppercase (keep)
             "1\trs001\t0\t2000\tc\tA"});  // c lowercase, A uppercase (reverse)

        auto bim_path = files.create_named_text_file("test.bim", bim_content);
        auto bed_path = bim_path;
        bed_path.replace_extension(".bed");

        SnpMatcher matcher(effects);
        auto match_plan = matcher.match(bed_path);

        REQUIRE(match_plan.size() == 2);
        REQUIRE(match_plan[0].type == MatchType::keep);     // a/A and c/C match
        REQUIRE(match_plan[1].type == MatchType::reverse);  // c/c and A/A match
    }
}

TEST_CASE("SnpMatcher - Edge cases", "[predictor][snp_matcher]")
{
    FileFixture files;

    SECTION("Happy path - single SNP in both files")
    {
        auto effects = create_snp_effects(
            files,
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs001\tA\tC\t0.25\t0.123\t0.045"});

        std::string bim_content
            = create_bim_content({"1\trs001\t0\t1000\tA\tC"});

        auto bim_path = files.create_named_text_file("test.bim", bim_content);
        auto bed_path = bim_path;
        bed_path.replace_extension(".bed");

        SnpMatcher matcher(effects);
        auto match_plan = matcher.match(bed_path);

        REQUIRE(match_plan.size() == 1);
        REQUIRE(match_plan[0].type == MatchType::keep);
        REQUIRE(match_plan[0].target_col == 0);
    }

    SECTION("Happy path - empty SNP effect file with non-empty BIM")
    {
        std::string empty_content = "ID\tA1\tA2\tA1Frq\tAdd\tDom\n";
        auto snp_effect_path
            = files.create_text_file(empty_content, ".snp.eff");
        gelex::detail::SnpEffectLoader loader(snp_effect_path);
        auto effects = std::move(loader).take_effects();

        std::string bim_content = create_bim_content(
            {"1\trs001\t0\t1000\tA\tC", "1\trs002\t0\t2000\tT\tG"});

        auto bim_path = files.create_named_text_file("test.bim", bim_content);
        auto bed_path = bim_path;
        bed_path.replace_extension(".bed");

        SnpMatcher matcher(effects);
        auto match_plan = matcher.match(bed_path);

        REQUIRE(match_plan.size() == 2);
        REQUIRE(match_plan[0].type == MatchType::skip);
        REQUIRE(match_plan[0].target_col == -1);
        REQUIRE(match_plan[1].type == MatchType::skip);
        REQUIRE(match_plan[1].target_col == -1);
    }

    SECTION("Happy path - SNP effect file without Dom column")
    {
        auto effects = create_snp_effects(
            files,
            "ID\tA1\tA2\tA1Frq\tAdd",
            {"rs001\tA\tC\t0.25\t0.123", "rs002\tT\tG\t0.75\t-0.456"});

        std::string bim_content = create_bim_content(
            {"1\trs001\t0\t1000\tA\tC", "1\trs002\t0\t2000\tT\tG"});

        auto bim_path = files.create_named_text_file("test.bim", bim_content);
        auto bed_path = bim_path;
        bed_path.replace_extension(".bed");

        SnpMatcher matcher(effects);
        auto match_plan = matcher.match(bed_path);

        REQUIRE(match_plan.size() == 2);
        REQUIRE(match_plan[0].type == MatchType::keep);
        REQUIRE(match_plan[0].target_col == 0);
        REQUIRE(match_plan[1].type == MatchType::keep);
        REQUIRE(match_plan[1].target_col == 1);
    }
}
