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

#include "snp_matcher.h"

#include <Eigen/Core>
#include <ranges>

#include "gelex/exception.h"

#include "../src/data/loader/bim_loader.h"
#include "gelex/types/snp_info.h"

namespace gelex
{

SnpMatcher::SnpMatcher(const SnpEffects& effects) : effects_(&effects) {}

MatchPlan SnpMatcher::match(const std::filesystem::path& predict_bed_path) const
{
    auto bim_path = predict_bed_path;
    bim_path.replace_extension(".bim");

    if (!std::filesystem::exists(bim_path))
    {
        throw FileNotFoundException(bim_path);
    }

    detail::BimLoader bim_loader(bim_path);

    SnpEffects predict_snp_meta = std::move(bim_loader).take_info();

    MatchPlan match_plan;
    match_plan.plan.resize(predict_snp_meta.size());
    match_plan.num_snp_in_effect = static_cast<Eigen::Index>(effects_->size());

    for (auto&& [meta, info] :
         std::views::zip(predict_snp_meta, match_plan.plan))
    {
        const auto snp_index = effects_->find_index(meta.id);

        if (snp_index)
        {
            const auto& model_meta = (*effects_)[*snp_index];
            info.type = determine_match_type(model_meta, meta);

            if (info.type != MatchType::skip)
            {
                info.target_col = *snp_index;
            }
        }
    }

    return match_plan;
}

constexpr char SnpMatcher::normalize_allele(char allele) noexcept
{
    return allele & 0xDF;
}

MatchType SnpMatcher::determine_match_type(
    const SnpMeta& model,
    const SnpMeta& predict) noexcept
{
    const char model_a1 = normalize_allele(model.A1);
    const char model_a2 = normalize_allele(model.A2);
    const char predict_a1 = normalize_allele(predict.A1);
    const char predict_a2 = normalize_allele(predict.A2);

    if (model_a1 == predict_a1 && model_a2 == predict_a2)
    {
        return MatchType::keep;
    }

    if (model_a1 == predict_a2 && model_a2 == predict_a1)
    {
        return MatchType::reverse;
    }

    return MatchType::skip;
}

}  // namespace gelex
