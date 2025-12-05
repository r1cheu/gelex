#ifndef GELEX_PREDICTOR_SNP_MATCHER_H_
#define GELEX_PREDICTOR_SNP_MATCHER_H_

#include <filesystem>
#include <vector>

#include <Eigen/Core>
#include "../src/data/loader/snp_effect_loader.h"

namespace gelex::detail
{
struct SnpInfo;
}
namespace gelex
{

enum class MatchType : uint8_t
{
    keep,
    reverse,
    skip
};

struct MatchInfo
{
    MatchType type = MatchType::skip;
    Eigen::Index target_col = -1;
};

using MatchPlan = std::vector<MatchInfo>;

class SnpMatcher
{
   public:
    explicit SnpMatcher(const std::filesystem::path& snp_effect_path);

    [[nodiscard]] MatchPlan match(
        const std::filesystem::path& predict_bed_path) const;

    SnpEffects take_snp_effects() && { return std::move(snp_effects_); }

   private:
    static constexpr char normalize_allele(char allele) noexcept;

    static MatchType determine_match_type(
        const SnpEffect& model,
        const detail::SnpInfo& predict) noexcept;

    SnpEffects snp_effects_;
};

}  // namespace gelex

#endif  // guard GELEX_PREDICTOR_SNP_MATCHER_H_
