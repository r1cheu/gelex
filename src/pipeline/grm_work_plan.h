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

#ifndef GELEX_PIPELINE_GRM_WORK_PLAN_H_
#define GELEX_PIPELINE_GRM_WORK_PLAN_H_

#include <filesystem>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "gelex/types/freq_effect.h"
#include "gelex/types/snp_info.h"

namespace gelex
{

struct GrmWorkItem
{
    std::vector<std::pair<Eigen::Index, Eigen::Index>> ranges;
    bool is_additive;
    std::string output_name;
};

class GrmNormalPlan
{
   public:
    GrmNormalPlan(const std::filesystem::path& bim_path, freq::GrmType mode);

    auto items() const -> const std::vector<GrmWorkItem>&;
    auto total_work() const -> Eigen::Index;
    auto output_pattern(std::string_view out_prefix) const -> std::string;

   private:
    std::vector<GrmWorkItem> items_;
    Eigen::Index total_work_ = 0;
    std::string task_pattern_;
};

class GrmLocoPlan
{
   public:
    GrmLocoPlan(const std::filesystem::path& bim_path, freq::GrmType mode);

    auto items() const -> const std::vector<GrmWorkItem>&;
    auto total_work() const -> Eigen::Index;
    auto output_pattern(std::string_view out_prefix) const -> std::string;

   private:
    struct ChrRange
    {
        std::string name;
        std::vector<std::pair<Eigen::Index, Eigen::Index>> ranges;
        Eigen::Index total_snps;
    };

    static auto build_loco_ranges(const SnpEffects& snp_effects)
        -> std::vector<ChrRange>;

    std::vector<GrmWorkItem> items_;
    Eigen::Index total_work_ = 0;
    size_t num_groups_ = 0;
    std::string task_pattern_;
};

}  // namespace gelex

#endif  // GELEX_PIPELINE_GRM_WORK_PLAN_H_
