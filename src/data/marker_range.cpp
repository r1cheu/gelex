// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

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
