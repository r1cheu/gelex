#ifndef GELEX_PREDICTOR_SNP_MATCHER_H_
#define GELEX_PREDICTOR_SNP_MATCHER_H_

#include <filesystem>
#include <vector>

#include <Eigen/Core>
#include "gelex/types/snp_info.h"

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

struct MatchPlan
{
    std::vector<MatchInfo> plan;
    Eigen::Index num_snp_in_effect = 0;

    MatchInfo& operator[](Eigen::Index idx) { return plan[idx]; }
    const MatchInfo& operator[](Eigen::Index idx) const { return plan[idx]; }
    void clear() noexcept
    {
        plan.clear();
        num_snp_in_effect = 0;
    }

    [[nodiscard]] size_t size() const noexcept { return plan.size(); }
};

class SnpMatcher
{
   public:
    explicit SnpMatcher(const SnpEffects& effects);

    [[nodiscard]] MatchPlan match(
        const std::filesystem::path& predict_bed_path) const;

   private:
    static constexpr char normalize_allele(char allele) noexcept;

    static MatchType determine_match_type(
        const SnpMeta& model,
        const SnpMeta& predict) noexcept;

    const SnpEffects* effects_;
};

}  // namespace gelex

#endif  // guard GELEX_PREDICTOR_SNP_MATCHER_H_
