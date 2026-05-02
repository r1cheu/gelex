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

#ifndef GELEX_PREDICT_SNP_ALIGNMENT_H_
#define GELEX_PREDICT_SNP_ALIGNMENT_H_

#include <cstddef>
#include <string>
#include <vector>

#include "gelex/data/dataframe/dataframe.h"

namespace gelex
{

inline constexpr double kMaxSnpMissingRatio = 0.2;

struct SnpAlignment
{
    std::vector<std::ptrdiff_t> column_map;  // -1 = missing or mismatched
    Eigen::Index num_missing{};
    Eigen::Index num_mismatched{};
};

[[nodiscard]] auto build_snp_alignment(
    const dataframe::DataFrame<std::string>& snp_effects,
    const dataframe::DataFrame<std::string>& bim_df) -> SnpAlignment;

}  // namespace gelex

#endif  // GELEX_PREDICT_SNP_ALIGNMENT_H_
