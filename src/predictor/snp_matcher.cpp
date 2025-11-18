#include "snp_matcher.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

namespace gelex
{

SnpMatcher::SnpMatcher(
    const std::unordered_map<std::string, SnpEffect>& model_snps)
    : model_snps_(model_snps)
{
    // Build reverse index for model SNPs
    int index = 0;
    for (const auto& [id, snp] : model_snps_)
    {
        model_indices_[id] = index++;
    }
}

std::expected<MatchingPlan, Error> SnpMatcher::match(
    const std::string& user_marker_info_file)
{
    auto user_snps_result = load_user_snps(user_marker_info_file);
    if (!user_snps_result.has_value())
    {
        return std::unexpected(user_snps_result.error());
    }

    return match(user_snps_result.value());
}

std::expected<MatchingPlan, Error> SnpMatcher::match(
    const std::vector<SnpEffect>& user_snps)
{
    MatchingPlan plan;
    plan.total_snps_in_user_file = static_cast<int>(user_snps.size());
    plan.total_snps_in_model = static_cast<int>(model_snps_.size());
    plan.matched_snps = 0;
    plan.flipped_snps = 0;
    plan.missing_snps = 0;

    for (int user_index = 0; user_index < static_cast<int>(user_snps.size());
         ++user_index)
    {
        const auto& user_snp = user_snps[user_index];

        // Look for matching SNP ID in model
        auto model_it = model_snps_.find(user_snp.id);
        if (model_it == model_snps_.end())
        {
            // SNP ID not found in model
            ++plan.missing_snps;
            continue;
        }

        const auto& model_snp = model_it->second;
        int model_index = model_indices_.at(user_snp.id);

        // Determine allele action
        AlleleAction action = determine_allele_action(model_snp, user_snp);

        if (action != AlleleAction::MISSING)
        {
            SnpMatch match;
            match.user_file_index = user_index;
            match.model_param_index = model_index;
            match.action = action;
            plan.matches.push_back(match);

            if (action == AlleleAction::FLIP)
            {
                ++plan.flipped_snps;
            }
            else
            {
                ++plan.matched_snps;
            }
        }
        else
        {
            ++plan.missing_snps;
        }
    }

    return plan;
}

AlleleAction SnpMatcher::determine_allele_action(
    const SnpEffect& model_snp,
    const SnpEffect& user_snp)
{
    // Convert alleles to uppercase for case-insensitive comparison
    char model_a1 = std::toupper(model_snp.a1);
    char model_a2 = std::toupper(model_snp.a2);
    char user_a1 = std::toupper(user_snp.a1);
    char user_a2 = std::toupper(user_snp.a2);

    // Check for direct match
    if (model_a1 == user_a1 && model_a2 == user_a2)
    {
        return AlleleAction::KEEP;
    }

    // Check for flipped match (A1/A2 swapped)
    if (model_a1 == user_a2 && model_a2 == user_a1)
    {
        return AlleleAction::FLIP;
    }

    // No valid match found
    return AlleleAction::MISSING;
}

bool SnpMatcher::are_complementary(char a1, char a2)
{
    // Convert to uppercase for case-insensitive comparison
    char upper_a1 = std::toupper(a1);
    char upper_a2 = std::toupper(a2);

    // Complementary base pairs
    return (upper_a1 == 'A' && upper_a2 == 'T')
           || (upper_a1 == 'T' && upper_a2 == 'A')
           || (upper_a1 == 'C' && upper_a2 == 'G')
           || (upper_a1 == 'G' && upper_a2 == 'C');
}

bool SnpMatcher::alleles_match(
    char model_a1,
    char model_a2,
    char user_a1,
    char user_a2)
{
    // Direct match
    if (model_a1 == user_a1 && model_a2 == user_a2)
    {
        return true;
    }

    // Flipped match
    if (model_a1 == user_a2 && model_a2 == user_a1)
    {
        return true;
    }

    // Complementary direct match
    if (are_complementary(model_a1, user_a1)
        && are_complementary(model_a2, user_a2))
    {
        return true;
    }

    // Complementary flipped match
    if (are_complementary(model_a1, user_a2)
        && are_complementary(model_a2, user_a1))
    {
        return true;
    }

    return false;
}

std::expected<std::vector<SnpEffect>, Error> SnpMatcher::load_user_snps(
    const std::string& marker_file)
{
    std::filesystem::path path(marker_file);
    std::string extension = path.extension().string();

    if (extension == ".bim")
    {
        return load_bim_file(marker_file);
    }
    else
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidFile,
                "Unsupported marker file format: " + extension
                    + ". Supported formats: .bim, .map"});
    }
}

std::expected<std::vector<SnpEffect>, Error> SnpMatcher::load_bim_file(
    const std::string& bim_file)
{
    std::ifstream file(bim_file);
    if (!file.is_open())
    {
        return std::unexpected(
            Error{
                ErrorCode::FileNotFound,
                "Failed to open BIM file: " + bim_file});
    }

    std::vector<SnpEffect> snps;
    std::string line;
    int line_number = 0;

    while (std::getline(file, line))
    {
        ++line_number;
        std::istringstream iss(line);
        std::string chrom, id, cm, pos, a1, a2;

        if (!(iss >> chrom >> id >> cm >> pos >> a1 >> a2))
        {
            return std::unexpected(
                Error{
                    ErrorCode::InvalidData,
                    "Failed to parse BIM file line "
                        + std::to_string(line_number) + ": " + line});
        }

        SnpEffect snp;
        snp.chrom = chrom;
        snp.id = id;
        snp.position = std::stoi(pos);

        // Convert allele strings to single characters
        // Take first character of allele string
        snp.a1 = a1.empty() ? '0' : a1[0];
        snp.a2 = a2.empty() ? '0' : a2[0];

        // Set default values for effect-related fields
        snp.p_freq = 0.0;
        snp.add_effect = 0.0;
        snp.dom_effect = 0.0;

        snps.push_back(snp);
    }

    return snps;
}
}  // namespace gelex
