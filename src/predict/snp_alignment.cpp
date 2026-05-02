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

#include "gelex/predict/snp_alignment.h"

#include <cstddef>
#include <string>

#include <Eigen/Core>

#include "gelex/data/dataframe/dataframe.h"

namespace gelex::predict
{

auto build_snp_alignment(
    const dataframe::DataFrame<std::string>& snp_effects,
    const dataframe::DataFrame<std::string>& bim_df) -> SnpAlignment
{
    const auto& eff_index = snp_effects.index();
    const auto& bim_index = bim_df.index();
    auto eff_a1 = snp_effects["A1"].as<std::string>();
    auto eff_a2 = snp_effects["A2"].as<std::string>();
    auto bim_a1 = bim_df["A1"].as<std::string>();
    auto bim_a2 = bim_df["A2"].as<std::string>();

    SnpAlignment result;
    result.column_map.reserve(eff_index.size());

    for (std::size_t i = 0; i < eff_index.size(); ++i)
    {
        const auto& key = eff_index.keys()[i];
        if (!bim_index.contains(key))
        {
            result.column_map.push_back(-1);
            result.num_missing++;
            continue;
        }
        auto bim_row = bim_index.at(key);
        if (eff_a1[i] == bim_a1[bim_row] && eff_a2[i] == bim_a2[bim_row])
        {
            result.column_map.push_back(static_cast<std::ptrdiff_t>(bim_row));
        }
        else
        {
            result.column_map.push_back(-1);
            result.num_mismatched++;
        }
    }
    return result;
}

}  // namespace gelex::predict
