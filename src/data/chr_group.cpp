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

#include "gelex/data/chr_group.h"
#include <Eigen/Core>
#include <cstddef>
#include <string>
#include <vector>
#include "gelex/data/dataframe/dataframe.h"

namespace gelex
{

auto build_chr_groups(
    bool do_loco,
    const dataframe::DataFrame<std::string>& bim) -> std::vector<ChrGroup>
{
    std::vector<ChrGroup> groups;
    auto num_snps = static_cast<Eigen::Index>(bim.rows());

    if (do_loco)
    {
        auto chrom = bim["chrom"].as<std::string>();
        std::string current_chr;
        Eigen::Index range_start = 0;

        for (Eigen::Index i = 0; i < num_snps; ++i)
        {
            if (chrom[static_cast<std::size_t>(i)] != current_chr)
            {
                if (!current_chr.empty())
                {
                    groups.push_back(
                        {current_chr, {{range_start, i}}, i - range_start});
                }
                current_chr = chrom[static_cast<std::size_t>(i)];
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
