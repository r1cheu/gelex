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

#include "gelex/pipeline/grm/grm_engine.h"

#include <filesystem>
#include <functional>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/data/genotype/genotype_method_dispatch.h"
#include "gelex/data/grm/grm.h"
#include "gelex/data/grm/grm_bin_writer.h"
#include "gelex/data/grm/grm_id_writer.h"
#include "gelex/data/loader/bim_loader.h"
#include "gelex/infra/logging/grm_event.h"
#include "gelex/infra/utils/utils.h"
#include "gelex/types/freq_effect.h"
#include "gelex/types/snp_info.h"

namespace gelex
{

namespace
{

struct ChrRange
{
    std::string name;
    std::vector<std::pair<Eigen::Index, Eigen::Index>> ranges;
    Eigen::Index total_snps;
};

struct GrmTask
{
    std::string name;
    bool is_additive;
};

template <typename Method>
auto dispatch_grm(
    GRM& grm,
    const std::vector<std::pair<Eigen::Index, Eigen::Index>>& ranges,
    int chunk_size,
    bool additive,
    const std::function<void(Eigen::Index, Eigen::Index)>& progress_callback
    = nullptr) -> GrmResult
{
    if (additive)
    {
        return grm.compute<typename Method::Additive>(
            ranges, chunk_size, progress_callback);
    }
    return grm.compute<typename Method::Dominant>(
        ranges, chunk_size, progress_callback);
}

auto compute_grm_with_method(
    GRM& grm,
    const std::vector<std::pair<Eigen::Index, Eigen::Index>>& ranges,
    GenotypeProcessMethod method,
    int chunk_size,
    bool additive,
    const std::function<void(Eigen::Index, Eigen::Index)>& progress_callback
    = nullptr) -> GrmResult
{
    auto visitor = [&]<typename MethodBundle>() -> GrmResult
    {
        return dispatch_grm<MethodBundle>(
            grm, ranges, chunk_size, additive, progress_callback);
    };
    return visit_genotype_method(method, visitor);
}

auto write_grm_files(
    const GrmResult& result,
    const std::vector<std::string>& sample_ids,
    const std::string& out_prefix) -> std::vector<std::string>
{
    auto bin_path = out_prefix + ".bin";
    auto id_path = out_prefix + ".id";

    GrmBinWriter(bin_path).write(result.grm);
    GrmIdWriter(id_path).write(sample_ids);

    return {bin_path, id_path};
}

auto build_chr_ranges(bool do_loco, const SnpEffects& snp_effects)
    -> std::vector<ChrRange>
{
    std::vector<ChrRange> groups;
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

auto build_tasks(gelex::freq::GrmType mode) -> std::vector<GrmTask>
{
    std::vector<GrmTask> tasks;
    if (mode != gelex::freq::GrmType::D)
    {
        tasks.push_back({"add", true});
    }
    if (mode != gelex::freq::GrmType::A)
    {
        tasks.push_back({"dom", false});
    }
    return tasks;
}

auto notify(const GrmObserver& observer, GrmEvent event) -> void
{
    if (observer)
    {
        observer(event);
    }
}

}  // namespace

GrmEngine::GrmEngine(Config config) : config_(std::move(config)) {}

auto GrmEngine::compute(const GrmObserver& observer) -> void
{
    GRM grm(config_.bed_path);
    const auto& sample_ids = grm.sample_ids();
    auto num_snps = static_cast<size_t>(grm.num_snps());

    notify(
        observer,
        GrmDataLoadedEvent{
            .num_samples = sample_ids.size(),
            .num_snps = num_snps,
        });

    auto bim_path = config_.bed_path;
    bim_path.replace_extension(".bim");
    detail::BimLoader bim_loader(bim_path);

    auto groups = build_chr_ranges(config_.do_loco, bim_loader.info());
    auto tasks = build_tasks(config_.mode);

    size_t total_work = 0;
    for (const auto& group : groups)
    {
        total_work += static_cast<size_t>(group.total_snps) * tasks.size();
    }

    SmoothEtaCalculator eta(total_work);
    std::vector<std::string> generated_files;
    size_t completed_base = 0;

    for (const auto& group : groups)
    {
        for (const auto& task : tasks)
        {
            auto progress_callback = [&](Eigen::Index current, Eigen::Index)
            {
                auto current_total
                    = completed_base + static_cast<size_t>(current);
                notify(
                    observer,
                    GrmProgressEvent{
                        .current = current_total,
                        .total = total_work,
                        .done = false,
                    });
            };

            auto result = compute_grm_with_method(
                grm,
                group.ranges,
                config_.method,
                config_.chunk_size,
                task.is_additive,
                progress_callback);

            auto suffix = config_.do_loco ? fmt::format(".chr{}", group.name)
                                          : std::string{};
            auto path
                = fmt::format("{}.{}{}", config_.out_prefix, task.name, suffix);

            auto files = write_grm_files(result, sample_ids, path);
            generated_files.insert(
                generated_files.end(), files.begin(), files.end());

            completed_base += static_cast<size_t>(group.total_snps);
        }
    }

    notify(
        observer,
        GrmProgressEvent{
            .current = total_work,
            .total = total_work,
            .done = true,
        });

    auto task_pattern
        = tasks.size() == 1 ? tasks[0].name : std::string("{add|dom}");
    std::string file_pattern;
    if (config_.do_loco)
    {
        file_pattern = fmt::format(
            "{}.{}.chr{{1..{}}}.{{bin|id}}",
            config_.out_prefix,
            task_pattern,
            groups.size());
    }
    else
    {
        file_pattern
            = fmt::format("{}.{}.{{bin|id}}", config_.out_prefix, task_pattern);
    }

    notify(
        observer,
        GrmFilesWrittenEvent{
            .file_paths = generated_files,
            .time_elapsed = eta.total_time_consumed(),
            .output_dir = std::filesystem::absolute(
                              std::filesystem::path(config_.out_prefix))
                              .parent_path()
                              .string(),
            .file_pattern = file_pattern,
        });
}

}  // namespace gelex
