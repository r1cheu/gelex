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

#include "gelex/pipeline/grm_engine.h"

#include <filesystem>
#include <string>
#include <vector>

#include <fmt/format.h>

#include "gelex/data/grm/grm.h"
#include "gelex/data/grm/grm_bin_writer.h"
#include "gelex/data/grm/grm_id_writer.h"
#include "gelex/infra/logging/grm_event.h"
#include "gelex/infra/utils/utils.h"
#include "pipeline/grm_work_plan.h"

namespace gelex
{

namespace
{

auto write_grm_files(
    const GrmResult& result,
    const std::vector<std::string>& sample_ids,
    const std::string& out_prefix) -> void
{
    GrmBinWriter(out_prefix + ".bin").write(result.grm);
    GrmIdWriter(out_prefix + ".id").write(sample_ids);
}

auto notify(const GrmObserver& observer, const GrmEvent& event) -> void
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

    notify(
        observer,
        GrmDataLoadedEvent{
            .num_samples = sample_ids.size(),
            .num_snps = static_cast<size_t>(grm.num_snps()),
        });

    auto dispatch_grm = [&](const auto& ranges, bool is_additive) -> GrmResult
    {
        if (is_additive)
        {
            return grm.compute<GeneticEffectType::Add>(
                config_.method, ranges, config_.chunk_size, observer);
        }
        return grm.compute<GeneticEffectType::Dom>(
            config_.method, ranges, config_.chunk_size, observer);
    };

    auto bim_path = config_.bed_path;
    bim_path.replace_extension(".bim");

    auto run_plan = [&](auto& plan)
    {
        const auto total_snps = static_cast<size_t>(plan.total_work());
        notify(observer, GrmComputeStartedEvent{.total_snps = total_snps});

        SmoothEtaCalculator eta(total_snps);

        for (const auto& item : plan.items())
        {
            auto result = dispatch_grm(item.ranges, item.is_additive);

            auto path
                = fmt::format("{}.{}", config_.out_prefix, item.output_name);
            write_grm_files(result, sample_ids, path);
        }

        notify(
            observer,
            GrmProgressEvent{
                .current = total_snps,
                .total = total_snps,
                .done = true,
            });

        notify(
            observer,
            GrmFilesWrittenEvent{
                .num_files = plan.items().size() * 2,
                .output_dir = std::filesystem::absolute(
                                  std::filesystem::path(config_.out_prefix))
                                  .parent_path()
                                  .string(),
                .file_pattern = plan.output_pattern(config_.out_prefix),
            });
    };

    if (config_.do_loco)
    {
        GrmLocoPlan plan(bim_path, config_.mode);
        run_plan(plan);
    }
    else
    {
        GrmNormalPlan plan(bim_path, config_.mode);
        run_plan(plan);
    }
}

}  // namespace gelex
