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

namespace gelex
{

auto build_snp_alignment(
    const DataFrame<std::string>& snp_effects,
    const DataFrame<std::string>& bim_df) -> AlignmentPlan
{
    const auto& eff_index = snp_effects.index();
    const auto& bim_index = bim_df.index();
    const auto eff_a1 = snp_effects["A1"].as<std::string>();
    const auto eff_a2 = snp_effects["A2"].as<std::string>();
    const auto bim_a1 = bim_df["A1"].as<std::string>();
    const auto bim_a2 = bim_df["A2"].as<std::string>();

    AlignmentPlan plan;
    plan.source_col.reserve(eff_index.size());
    plan.train_pos.reserve(eff_index.size());
    plan.flip.reserve(eff_index.size());

    for (std::size_t i = 0; i < eff_index.size(); ++i)
    {
        const auto& key = eff_index.keys()[i];
        const auto* row = bim_index.find(key);
        if (row == nullptr)
        {
            plan.missing_pos.push_back(static_cast<Eigen::Index>(i));
            ++plan.num_absent;
            continue;
        }

        if (eff_a1[i] == bim_a1[*row] && eff_a2[i] == bim_a2[*row])
        {
            plan.source_col.push_back(static_cast<Eigen::Index>(*row));
            plan.train_pos.push_back(static_cast<Eigen::Index>(i));
            plan.flip.push_back(0);
            ++plan.num_same;
        }
        else if (eff_a1[i] == bim_a2[*row] && eff_a2[i] == bim_a1[*row])
        {
            plan.source_col.push_back(static_cast<Eigen::Index>(*row));
            plan.train_pos.push_back(static_cast<Eigen::Index>(i));
            plan.flip.push_back(1);
            ++plan.num_flip;
        }
        else
        {
            plan.missing_pos.push_back(static_cast<Eigen::Index>(i));
            ++plan.num_incompatible;
        }
    }

    return plan;
}

}  // namespace gelex
