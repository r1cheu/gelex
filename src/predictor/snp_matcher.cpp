#include "snp_matcher.h"

#include <cctype>
#include <filesystem>
#include <span>
#include <unordered_map>

#include "data/loader.h"
#include "predictor/snp_effect_processor.h"

namespace gelex
{

auto SnpMatcher::create(const std::filesystem::path& bed_prefix)
    -> std::expected<SnpMatcher, Error>
{
    std::filesystem::path bim_path = bed_prefix;
    bim_path.replace_extension(".bim");

    auto bim_loader_result = detail::BimLoader::create(bim_path);
    if (!bim_loader_result.has_value())
    {
        return std::unexpected(bim_loader_result.error());
    }
    SnpMatcher matcher;
    matcher.predict_snp_meta_ = std::move(*bim_loader_result).take_meta();
    return matcher;
}

auto SnpMatcher::match(std::span<const SnpEffect> snp_effects) -> MatchInfo
{
    MatchInfo match_info;
    match_info.read_plan.resize(predict_snp_meta_.size(), MatchType::skip);

    // Build model column mapping
    for (size_t i = 0; i < snp_effects.size(); ++i)
    {
        match_info.model_column_map[snp_effects[i].meta.id] = i;
    }

    // Match prediction SNPs to model SNPs
    for (size_t i = 0; i < predict_snp_meta_.size(); ++i)
    {
        const auto& predict_snp = predict_snp_meta_[i];
        // Check if prediction SNP exists in model
        auto it = match_info.model_column_map.find(predict_snp.id);
        if (it != match_info.model_column_map.end())
        {
            const auto& model_snp = snp_effects[it->second];
            match_info.read_plan[i]
                = determine_match_type(model_snp.meta, predict_snp);

            // Add to prediction column order
            match_info.prediction_column_order.push_back(it->second);
        }
    }

    return match_info;
}

char SnpMatcher::normalize_allele(char allele)
{
    return static_cast<char>(std::toupper(allele));
}

MatchType SnpMatcher::determine_match_type(
    const detail::SnpMeta& model,
    const detail::SnpMeta& predict)
{
    // Normalize alleles for case-insensitive comparison
    char model_a1 = normalize_allele(model.a1);
    char model_a2 = normalize_allele(model.a2);
    char predict_a1 = normalize_allele(predict.a1);
    char predict_a2 = normalize_allele(predict.a2);

    // Direct match: alleles match exactly
    if (model_a1 == predict_a1 && model_a2 == predict_a2)
    {
        return MatchType::keep;
    }

    // Reverse match: alleles are swapped
    if (model_a1 == predict_a2 && model_a2 == predict_a1)
    {
        return MatchType::reverse;
    }

    // No match found
    return MatchType::skip;
}

}  // namespace gelex
