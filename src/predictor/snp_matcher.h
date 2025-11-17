#pragma once

#include <expected>
#include <string>
#include <unordered_map>
#include <vector>

#include "gelex/error.h"

namespace gelex
{

/// Unified SNP information structure combining metadata and effect data
struct SnpInfo
{
    std::string id;
    std::string chrom;
    int position;
    char a1;
    char a2;
    double p_freq;      // A1 allele frequency
    double add_effect;  // Additive effect coefficient
    double dom_effect;  // Dominant effect coefficient
};

/// Action to take for allele matching
enum class AlleleAction : uint8_t
{
    KEEP,    // Keep genotype as-is
    FLIP,    // Flip genotype (0->2, 2->0)
    MISSING  // Mark as missing (no valid match)
};

/// Individual SNP match information
struct SnpMatch
{
    int user_file_index;    // User file SNP index
    int model_param_index;  // Model parameter SNP index
    AlleleAction action;    // Action to take for allele matching
};

/// Complete matching plan between user SNPs and model SNPs
struct MatchingPlan
{
    std::vector<SnpMatch> matches;
    int total_snps_in_user_file;
    int total_snps_in_model;
    int matched_snps;
    int flipped_snps;
    int missing_snps;
};

/// SNP matcher for creating matching plans between model SNPs and user genotype
/// data Handles allele flips and strand issues without loading full genotype
/// matrices
class SnpMatcher
{
   public:
    /// Constructor with model SNPs
    /// @param model_snps Map of model SNP ID -> SnpInfo
    explicit SnpMatcher(
        const std::unordered_map<std::string, SnpInfo>& model_snps);

    /// Create matching plan between model SNPs and user marker file
    /// @param user_marker_info_file Path to user marker file (.bim, .map, or
    /// VCF header)
    /// @return Expected containing MatchingPlan, or Error on failure
    std::expected<MatchingPlan, Error> match(
        const std::string& user_marker_info_file);

    /// Create matching plan from vector of user SNP information
    /// @param user_snps Vector of user SNP information
    /// @return Expected containing MatchingPlan, or Error on failure
    std::expected<MatchingPlan, Error> match(
        const std::vector<SnpInfo>& user_snps);

    /// Check if alleles are complementary (A/T, C/G)
    /// @param a1 First allele
    /// @param a2 Second allele
    /// @return True if alleles are complementary
    static bool are_complementary(char a1, char a2);

    /// Check if alleles match directly or with flip
    /// @param model_a1 Model A1 allele
    /// @param model_a2 Model A2 allele
    /// @param user_a1 User A1 allele
    /// @param user_a2 User A2 allele
    /// @return True if alleles match (directly or with flip)
    static bool
    alleles_match(char model_a1, char model_a2, char user_a1, char user_a2);

   private:
    /// Model SNPs stored as ID -> SnpInfo map for fast lookup
    std::unordered_map<std::string, SnpInfo> model_snps_;

    /// Model SNP indices for reverse lookup
    std::unordered_map<std::string, int> model_indices_;

    /// Determine allele action for SNP matching
    /// @param model_snp Model SNP information
    /// @param user_snp User SNP information
    /// @return AlleleAction indicating what action to take
    AlleleAction determine_allele_action(
        const SnpInfo& model_snp,
        const SnpInfo& user_snp);

    /// Load user SNP information from marker file
    /// @param marker_file Path to marker file
    /// @return Expected containing vector of SnpInfo, or Error on failure
    std::expected<std::vector<SnpInfo>, Error> load_user_snps(
        const std::string& marker_file);

    /// Load user SNPs from BIM file
    /// @param bim_file Path to BIM file
    /// @return Expected containing vector of SnpInfo, or Error on failure
    std::expected<std::vector<SnpInfo>, Error> load_bim_file(
        const std::string& bim_file);
};

}  // namespace gelex
