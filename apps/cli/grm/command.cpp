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

#include "command.h"

#include <cstddef>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <CLI/CLI.hpp>
#include <Eigen/Core>

#include "cli/cli_helper.h"
#include "gelex/data/grm/grm.h"
#include "gelex/data/reader.h"
#include "gelex/data/writer.h"
#include "gelex/types/genetic_effect_type.h"
#include "reporter.h"

namespace
{

struct GrmTask
{
    std::string name;
    bool is_additive;
};

struct GrmWorkItem
{
    std::vector<std::pair<Eigen::Index, Eigen::Index>> ranges;
    bool is_additive;
    std::string output_name;
};

}  // namespace

auto grm_execute(CLI::App& cmd) -> int
{
    cli::GrmReporter reporter;

    auto threads = cmd.get_option("--threads")->as<int>();
    cli::setup_parallelization(threads);

    std::vector<gelex::GeneticMode> requested_effects;
    if (cmd.get_option("--add")->count() > 0)
    {
        requested_effects.push_back(gelex::GeneticMode::A);
    }
    if (cmd.get_option("--dom")->count() > 0)
    {
        requested_effects.push_back(gelex::GeneticMode::D);
    }
    if (requested_effects.empty())
    {
        requested_effects.push_back(gelex::GeneticMode::A);
    }

    auto method = cli::parse_genotype_method(
        cmd.get_option("--geno-method")->as<std::string>());
    auto chunk_size = cmd.get_option("--chunk-size")->as<int>();

    cli::GrmReporter::show_banner();
    cli::GrmReporter::show_config(
        fmt::format("{}", method),
        requested_effects,
        cmd.get_option("--loco")->count() > 0);

    gelex::GRM grm(cmd.get_option("--bfile")->as<std::string>());
    const auto& sample_ids = grm.sample_ids();
    const auto observer = reporter.as_observer();

    cli::GrmReporter::show_data_loaded(
        sample_ids.size(), static_cast<size_t>(grm.num_snps()));

    std::vector<GrmTask> tasks;
    for (auto effect : requested_effects)
    {
        if (effect == gelex::GeneticMode::A)
        {
            tasks.push_back({.name = "add", .is_additive = true});
        }
        if (effect == gelex::GeneticMode::D)
        {
            tasks.push_back({.name = "dom", .is_additive = false});
        }
    }

    std::vector<GrmWorkItem> items;
    Eigen::Index total_work = 0;
    const auto bfile_prefix = cmd.get_option("--bfile")->as<std::string>();
    auto bim = gelex::read_bim(bfile_prefix + ".bim");
    const auto num_snps = static_cast<Eigen::Index>(bim.rows());

    std::string task_pattern
        = tasks.size() == 1 ? tasks[0].name : std::string("{add|dom}");

    if (cmd.get_option("--loco")->count() > 0)
    {
        struct ChrRange
        {
            std::string name;
            std::vector<std::pair<Eigen::Index, Eigen::Index>> ranges;
            Eigen::Index total_snps;
        };

        std::vector<ChrRange> groups;
        auto chrom = bim["chrom"].as<std::string>();
        std::string current_chr;
        Eigen::Index range_start = 0;

        for (Eigen::Index i = 0; i < num_snps; ++i)
        {
            if (chrom[static_cast<std::size_t>(i)] != current_chr)
            {
                if (!current_chr.empty())
                {
                    groups.push_back(
                        {current_chr, {{range_start, i}}, i - range_start});
                }
                current_chr = chrom[static_cast<std::size_t>(i)];
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

        items.reserve(groups.size() * tasks.size());
        for (const auto& group : groups)
        {
            for (const auto& task : tasks)
            {
                items.push_back(
                    {.ranges = group.ranges,
                     .is_additive = task.is_additive,
                     .output_name
                     = fmt::format("{}.chr{}", task.name, group.name)});
                total_work += group.total_snps;
            }
        }
    }
    else
    {
        items.reserve(tasks.size());
        for (const auto& task : tasks)
        {
            items.push_back(
                {.ranges = {{0, num_snps}},
                 .is_additive = task.is_additive,
                 .output_name = task.name});
            total_work += num_snps;
        }
    }

    reporter.start_compute(static_cast<size_t>(total_work));

    for (const auto& item : items)
    {
        gelex::GrmResult result
            = item.is_additive ? grm.compute<gelex::GeneticMode::A>(
                                     method, item.ranges, chunk_size, observer)
                               : grm.compute<gelex::GeneticMode::D>(
                                     method, item.ranges, chunk_size, observer);

        gelex::write_grm(
            fmt::format(
                "{}.{}",
                cmd.get_option("--out")->as<std::string>(),
                item.output_name),
            result.grm,
            sample_ids);
    }

    reporter.finish_progress();

    auto output_pattern = cmd.get_option("--loco")->count() > 0
                              ? fmt::format(
                                    "{}.{}.chr{{1..{}}}.{{bin|id}}",
                                    cmd.get_option("--out")->as<std::string>(),
                                    task_pattern,
                                    items.size() / tasks.size())
                              : fmt::format(
                                    "{}.{}.{{bin|id}}",
                                    cmd.get_option("--out")->as<std::string>(),
                                    task_pattern);

    cli::GrmReporter::show_files_written(
        items.size() * 2,
        std::filesystem::absolute(
            std::filesystem::path(cmd.get_option("--out")->as<std::string>()))
            .parent_path()
            .string(),
        output_pattern);

    return 0;
}
