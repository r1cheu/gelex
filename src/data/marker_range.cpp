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

#include "gelex/data/marker_range.h"

#include <Eigen/Core>
#include <cstddef>
#include <string>
#include <vector>

#include "gelex/data/dataframe/dataframe.h"

namespace gelex
{

auto chromosome_ranges(const DataFrame<std::string>& bim)
    -> std::vector<MarkerRange>
{
    const auto num_snps = static_cast<Eigen::Index>(bim.rows());
    auto chrom = bim["CHR"].as<std::string>();

    std::vector<MarkerRange> ranges;
    std::string current;
    Eigen::Index range_start = 0;
    for (Eigen::Index i = 0; i < num_snps; ++i)
    {
        if (chrom[static_cast<std::size_t>(i)] != current)
        {
            if (!current.empty())
            {
                ranges.push_back({current, range_start, i});
            }
            current = chrom[static_cast<std::size_t>(i)];
            range_start = i;
        }
    }
    if (!current.empty())
    {
        ranges.push_back({current, range_start, num_snps});
    }
    return ranges;
}

}  // namespace gelex
