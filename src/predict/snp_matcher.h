/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef GELEX_PREDICT_SNP_MATCHER_H_
#define GELEX_PREDICT_SNP_MATCHER_H_

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

#endif  // guard GELEX_PREDICT_SNP_MATCHER_H_
