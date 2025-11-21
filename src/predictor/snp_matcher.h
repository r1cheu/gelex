#pragma once

#include <expected>
#include <filesystem>
#include <unordered_map>
#include <vector>

#include "data/loader.h"
#include "gelex/error.h"

namespace gelex
{
struct SnpEffect;

enum class MatchType : uint8_t
{
    keep,
    reverse,
    skip
};

/// Structure to track matching information between model and prediction SNPs
struct MatchInfo
{
    std::vector<MatchType> read_plan;  ///< Match type for each prediction SNP
    std::unordered_map<std::string, size_t>
        model_column_map;  ///< Model SNP ID to column index
    std::vector<size_t>
        prediction_column_order;  ///< Column order for prediction genotypes
};

class SnpMatcher
{
   public:
    static auto create(const std::filesystem::path& bed_prefix)
        -> std::expected<SnpMatcher, Error>;

    auto match(std::span<const SnpEffect> snp_effects) -> MatchInfo;

   private:
    SnpMatcher() = default;

    // Helper functions for allele matching
    static char normalize_allele(char allele);
    static MatchType determine_match_type(
        const detail::SnpMeta& model,
        const detail::SnpMeta& predict);

    std::vector<detail::SnpMeta> predict_snp_meta_;
};

}  // namespace gelex
