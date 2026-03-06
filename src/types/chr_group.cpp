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

#include "gelex/types/chr_group.h"

#include "gelex/types/snp_info.h"

namespace gelex
{

auto build_chr_groups(bool do_loco, const SnpEffects& snp_effects)
    -> std::vector<ChrGroup>
{
    std::vector<ChrGroup> groups;
    auto num_snps = static_cast<Eigen::Index>(snp_effects.size());

    if (do_loco)
    {
        std::string current_chr;
        Eigen::Index range_start = 0;

        for (Eigen::Index i = 0; i < num_snps; ++i)
        {
            if (snp_effects[i].chrom != current_chr)
            {
                if (!current_chr.empty())
                {
                    groups.push_back(
                        {current_chr, {{range_start, i}}, i - range_start});
                }
                current_chr = snp_effects[i].chrom;
                range_start = i;
            }
        }
        if (!current_chr.empty())
        {
            groups.push_back(
                {current_chr,
                 {{range_start, num_snps}},
                 num_snps - range_start});
        }
    }
    else
    {
        groups.push_back({"all", {{0, num_snps}}, num_snps});
    }
    return groups;
}

}  // namespace gelex
