#include <catch2/catch_test_macros.hpp>

#include "../src/predictor/snp_matcher.h"

using namespace gelex;

TEST_CASE("SnpMatcher - Basic SNP matching", "[snp_matcher]")
{
    // Create model SNPs
    std::unordered_map<std::string, SnpInfo> model_snps;

    model_snps["rs1"] = {"rs1", "1", 1000, 'A', 'T', 0.3, 0.1, 0.02};
    model_snps["rs2"] = {"rs2", "1", 2000, 'C', 'G', 0.4, 0.2, 0.03};
    model_snps["rs3"] = {"rs3", "1", 3000, 'A', 'C', 0.5, 0.3, 0.04};

    SnpMatcher matcher(model_snps);

    // Create user SNPs
    std::vector<SnpInfo> user_snps;
    user_snps.push_back({"rs1", "1", 1000, 'A', 'T', 0.0, 0.0, 0.0});  // Direct match
    user_snps.push_back({"rs2", "1", 2000, 'G', 'C', 0.0, 0.0, 0.0});  // Flipped match
    user_snps.push_back({"rs4", "1", 4000, 'A', 'G', 0.0, 0.0, 0.0});  // Missing SNP

    auto result = matcher.match(user_snps);
    REQUIRE(result.has_value());

    const auto& plan = result.value();

    REQUIRE(plan.total_snps_in_user_file == 3);
    REQUIRE(plan.total_snps_in_model == 3);
    REQUIRE(plan.matched_snps == 1);  // Direct match
    REQUIRE(plan.flipped_snps == 1);  // Flipped match
    REQUIRE(plan.missing_snps == 1);  // Missing SNP
    REQUIRE(plan.matches.size() == 2);

    // Check matches exist with correct actions
    REQUIRE(plan.matches.size() == 2);

    // Find the direct match
    auto direct_match = std::find_if(plan.matches.begin(), plan.matches.end(),
        [](const SnpMatch& match) { return match.action == AlleleAction::KEEP; });
    REQUIRE(direct_match != plan.matches.end());
    REQUIRE(direct_match->user_file_index == 0);

    // Find the flipped match
    auto flipped_match = std::find_if(plan.matches.begin(), plan.matches.end(),
        [](const SnpMatch& match) { return match.action == AlleleAction::FLIP; });
    REQUIRE(flipped_match != plan.matches.end());
    REQUIRE(flipped_match->user_file_index == 1);
}

TEST_CASE("SnpMatcher - Complementary allele matching", "[snp_matcher]")
{
    std::unordered_map<std::string, SnpInfo> model_snps;
    model_snps["rs1"] = {"rs1", "1", 1000, 'A', 'T', 0.3, 0.1, 0.02};
    model_snps["rs2"] = {"rs2", "1", 2000, 'C', 'G', 0.4, 0.2, 0.03};

    SnpMatcher matcher(model_snps);

    std::vector<SnpInfo> user_snps;
    user_snps.push_back({"rs1", "1", 1000, 'T', 'A', 0.0, 0.0, 0.0});  // Complementary direct
    user_snps.push_back({"rs2", "1", 2000, 'G', 'C', 0.0, 0.0, 0.0});  // Complementary flipped

    auto result = matcher.match(user_snps);
    REQUIRE(result.has_value());

    const auto& plan = result.value();

    REQUIRE(plan.matches.size() == 2);
    REQUIRE(plan.matched_snps == 0);  // No direct matches
    REQUIRE(plan.flipped_snps == 2);  // Both complementary cases are treated as flipped
}

TEST_CASE("SnpMatcher - Allele mismatch handling", "[snp_matcher]")
{
    std::unordered_map<std::string, SnpInfo> model_snps;
    model_snps["rs1"] = {"rs1", "1", 1000, 'A', 'T', 0.3, 0.1, 0.02};

    SnpMatcher matcher(model_snps);

    std::vector<SnpInfo> user_snps;
    user_snps.push_back({"rs1", "1", 1000, 'A', 'C', 0.0, 0.0, 0.0});  // Mismatched alleles

    auto result = matcher.match(user_snps);
    REQUIRE(result.has_value());

    const auto& plan = result.value();

    REQUIRE(plan.matches.empty());
    REQUIRE(plan.missing_snps == 1);
}

TEST_CASE("SnpMatcher - Empty model SNPs", "[snp_matcher]")
{
    std::unordered_map<std::string, SnpInfo> model_snps;
    SnpMatcher matcher(model_snps);

    std::vector<SnpInfo> user_snps;
    user_snps.push_back({"rs1", "1", 1000, 'A', 'T', 0.0, 0.0, 0.0});

    auto result = matcher.match(user_snps);
    REQUIRE(result.has_value());

    const auto& plan = result.value();

    REQUIRE(plan.total_snps_in_user_file == 1);
    REQUIRE(plan.total_snps_in_model == 0);
    REQUIRE(plan.matched_snps == 0);
    REQUIRE(plan.flipped_snps == 0);
    REQUIRE(plan.missing_snps == 1);
    REQUIRE(plan.matches.empty());
}

TEST_CASE("SnpMatcher - Empty user SNPs", "[snp_matcher]")
{
    std::unordered_map<std::string, SnpInfo> model_snps;
    model_snps["rs1"] = {"rs1", "1", 1000, 'A', 'T', 0.3, 0.1, 0.02};

    SnpMatcher matcher(model_snps);

    std::vector<SnpInfo> user_snps;

    auto result = matcher.match(user_snps);
    REQUIRE(result.has_value());

    const auto& plan = result.value();

    REQUIRE(plan.total_snps_in_user_file == 0);
    REQUIRE(plan.total_snps_in_model == 1);
    REQUIRE(plan.matched_snps == 0);
    REQUIRE(plan.flipped_snps == 0);
    REQUIRE(plan.missing_snps == 0);
    REQUIRE(plan.matches.empty());
}

TEST_CASE("SnpMatcher - Case insensitive allele matching", "[snp_matcher]")
{
    std::unordered_map<std::string, SnpInfo> model_snps;
    model_snps["rs1"] = {"rs1", "1", 1000, 'A', 'T', 0.3, 0.1, 0.02};

    SnpMatcher matcher(model_snps);

    std::vector<SnpInfo> user_snps;
    user_snps.push_back({"rs1", "1", 1000, 'a', 't', 0.0, 0.0, 0.0});  // Lowercase alleles

    auto result = matcher.match(user_snps);
    REQUIRE(result.has_value());

    const auto& plan = result.value();

    REQUIRE(plan.matches.size() == 1);
    REQUIRE(plan.matched_snps == 1);
    REQUIRE(plan.matches[0].action == AlleleAction::KEEP);
}

TEST_CASE("SnpMatcher - Mixed allele scenarios", "[snp_matcher]")
{
    std::unordered_map<std::string, SnpInfo> model_snps;
    model_snps["rs1"] = {"rs1", "1", 1000, 'A', 'T', 0.3, 0.1, 0.02};
    model_snps["rs2"] = {"rs2", "1", 2000, 'C', 'G', 0.4, 0.2, 0.03};
    model_snps["rs3"] = {"rs3", "1", 3000, 'A', 'C', 0.5, 0.3, 0.04};

    SnpMatcher matcher(model_snps);

    std::vector<SnpInfo> user_snps;
    user_snps.push_back({"rs1", "1", 1000, 'A', 'T', 0.0, 0.0, 0.0});  // Direct match
    user_snps.push_back({"rs1", "1", 1000, 'T', 'A', 0.0, 0.0, 0.0});  // Flipped match
    user_snps.push_back({"rs2", "1", 2000, 'G', 'C', 0.0, 0.0, 0.0});  // Complementary flipped
    user_snps.push_back({"rs3", "1", 3000, 'C', 'A', 0.0, 0.0, 0.0});  // Regular flipped
    user_snps.push_back({"rs4", "1", 4000, 'A', 'G', 0.0, 0.0, 0.0});  // Missing SNP

    auto result = matcher.match(user_snps);
    REQUIRE(result.has_value());

    const auto& plan = result.value();

    REQUIRE(plan.total_snps_in_user_file == 5);
    REQUIRE(plan.total_snps_in_model == 3);
    REQUIRE(plan.matched_snps == 1);  // Direct match only
    REQUIRE(plan.flipped_snps == 3);  // All flipped matches
    REQUIRE(plan.missing_snps == 1);  // Missing SNP
    REQUIRE(plan.matches.size() == 4);
}

TEST_CASE("SnpMatcher - are_complementary function", "[snp_matcher]")
{
    // Test complementary pairs
    REQUIRE(SnpMatcher::are_complementary('A', 'T') == true);
    REQUIRE(SnpMatcher::are_complementary('T', 'A') == true);
    REQUIRE(SnpMatcher::are_complementary('C', 'G') == true);
    REQUIRE(SnpMatcher::are_complementary('G', 'C') == true);

    // Test case insensitive
    REQUIRE(SnpMatcher::are_complementary('a', 't') == true);
    REQUIRE(SnpMatcher::are_complementary('A', 't') == true);

    // Test non-complementary pairs
    REQUIRE(SnpMatcher::are_complementary('A', 'A') == false);
    REQUIRE(SnpMatcher::are_complementary('A', 'C') == false);
    REQUIRE(SnpMatcher::are_complementary('C', 'T') == false);
    REQUIRE(SnpMatcher::are_complementary('G', 'A') == false);
}

TEST_CASE("SnpMatcher - alleles_match function", "[snp_matcher]")
{
    // Direct match
    REQUIRE(SnpMatcher::alleles_match('A', 'T', 'A', 'T') == true);

    // Flipped match
    REQUIRE(SnpMatcher::alleles_match('A', 'T', 'T', 'A') == true);

    // Complementary direct match
    REQUIRE(SnpMatcher::alleles_match('A', 'T', 'T', 'A') == true);

    // Complementary flipped match
    REQUIRE(SnpMatcher::alleles_match('A', 'T', 'A', 'T') == true);

    // Mismatched alleles
    REQUIRE(SnpMatcher::alleles_match('A', 'T', 'A', 'C') == false);
    REQUIRE(SnpMatcher::alleles_match('A', 'T', 'C', 'G') == false);
}