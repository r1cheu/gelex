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

#include "pipeline/grm_work_plan.h"

#include <string>
#include <vector>

#include <fmt/format.h>

#include "gelex/data/loader/bim_loader.h"

namespace gelex
{

namespace
{

struct GrmTask
{
    std::string name;
    bool is_additive;
};

auto build_tasks(freq::GrmType mode) -> std::vector<GrmTask>
{
    std::vector<GrmTask> tasks;
    if (mode != freq::GrmType::D)
    {
        tasks.push_back({"add", true});
    }
    if (mode != freq::GrmType::A)
    {
        tasks.push_back({"dom", false});
    }
    return tasks;
}

}  // namespace

auto GrmLocoPlan::build_loco_ranges(const SnpEffects& snp_effects)
    -> std::vector<ChrRange>
{
    std::vector<ChrRange> groups;
    auto num_snps = static_cast<Eigen::Index>(snp_effects.size());
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
            {current_chr, {{range_start, num_snps}}, num_snps - range_start});
    }
    return groups;
}

GrmNormalPlan::GrmNormalPlan(
    const std::filesystem::path& bim_path,
    freq::GrmType mode)
{
    detail::BimLoader bim_loader(bim_path);
    auto num_snps = static_cast<Eigen::Index>(bim_loader.info().size());
    auto tasks = build_tasks(mode);

    task_pattern_
        = tasks.size() == 1 ? tasks[0].name : std::string("{add|dom}");

    items_.reserve(tasks.size());
    for (const auto& task : tasks)
    {
        items_.push_back({
            .ranges = {{0, num_snps}},
            .is_additive = task.is_additive,
            .output_name = task.name,
        });
        total_work_ += num_snps;
    }
}

auto GrmNormalPlan::items() const -> const std::vector<GrmWorkItem>&
{
    return items_;
}

auto GrmNormalPlan::total_work() const -> Eigen::Index
{
    return total_work_;
}

auto GrmNormalPlan::output_pattern(std::string_view out_prefix) const
    -> std::string
{
    return fmt::format("{}.{}.{{bin|id}}", out_prefix, task_pattern_);
}

GrmLocoPlan::GrmLocoPlan(
    const std::filesystem::path& bim_path,
    freq::GrmType mode)
{
    detail::BimLoader bim_loader(bim_path);
    auto groups = build_loco_ranges(bim_loader.info());
    auto tasks = build_tasks(mode);

    num_groups_ = groups.size();
    task_pattern_
        = tasks.size() == 1 ? tasks[0].name : std::string("{add|dom}");

    items_.reserve(groups.size() * tasks.size());
    for (const auto& group : groups)
    {
        for (const auto& task : tasks)
        {
            items_.push_back({
                .ranges = group.ranges,
                .is_additive = task.is_additive,
                .output_name = fmt::format("{}.chr{}", task.name, group.name),
            });
            total_work_ += group.total_snps;
        }
    }
}

auto GrmLocoPlan::items() const -> const std::vector<GrmWorkItem>&
{
    return items_;
}

auto GrmLocoPlan::total_work() const -> Eigen::Index
{
    return total_work_;
}

auto GrmLocoPlan::output_pattern(std::string_view out_prefix) const
    -> std::string
{
    return fmt::format(
        "{}.{}.chr{{1..{}}}.{{bin|id}}",
        out_prefix,
        task_pattern_,
        num_groups_);
}

}  // namespace gelex
