// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_DATA_MARKER_RANGE_H_
#define GELEX_DATA_MARKER_RANGE_H_

#include <Eigen/Core>
#include <string>
#include <vector>

#include "gelex/data/dataframe/dataframe.h"

namespace gelex
{

struct MarkerRange
{
    std::string label;  // empty for whole-genome; chromosome name for per-chr
    Eigen::Index start;
    Eigen::Index end;
};

// Contiguous marker runs per chromosome, in bim order.
auto chromosome_ranges(const DataFrame<std::string>& bim)
    -> std::vector<MarkerRange>;

}  // namespace gelex

#endif  // GELEX_DATA_MARKER_RANGE_H_
