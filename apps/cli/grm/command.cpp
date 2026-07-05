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
#include <string_view>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "cli/cli_helper.h"
#include "gelex/data/bed.h"
#include "gelex/data/grm/grm.h"
#include "gelex/data/marker_range.h"
#include "gelex/data/reader.h"
#include "gelex/data/writer.h"
#include "gelex/types/genetic_mode.h"
#include "reporter.h"

namespace
{

auto mode_tag(gelex::GeneticMode mode) -> std::string_view
{
    return mode == gelex::GeneticMode::A ? "add" : "dom";
}

}  // namespace

auto grm_execute(const cli::GrmConfig& config) -> int
{
    cli::GrmReporter reporter;

    cli::setup_parallelization(config.threads);

    const auto modes = config.mode;
    const auto method = config.geno_method;
    const auto chunk_size = config.chunk_size;

    auto bed = gelex::open_bed(config.bfile);
    const auto sample_ids = bed.sample_index().keys();
    const auto observer = reporter.as_observer();

    cli::GrmReporter::show_data_loaded(
        sample_ids.size(), static_cast<size_t>(bed.num_snps()));

    std::vector<gelex::MarkerRange> ranges;
    if (config.loco)
    {
        auto bim = gelex::read_bim(config.bfile + ".bim");
        ranges = gelex::chromosome_ranges(bim);
    }
    else
    {
        ranges.push_back({std::string{}, 0, bed.num_snps()});
    }

    gelex::GrmBuilder builder(bed, modes, method, chunk_size, observer);
    builder.build(
        ranges,
        [&](const gelex::GrmMatrix& matrix)
        {
            const std::string name = matrix.label.empty()
                                         ? std::string{mode_tag(matrix.mode)}
                                         : fmt::format(
                                               "{}.chr{:02d}",
                                               mode_tag(matrix.mode),
                                               std::stoi(matrix.label));
            gelex::write_grm(
                fmt::format("{}.{}", config.out, name), matrix.grm, sample_ids);
        });

    reporter.finish_progress();

    const std::string task_pattern
        = modes.size() == 1 ? std::string{mode_tag(
                                  modes.contains(gelex::GeneticMode::A)
                                      ? gelex::GeneticMode::A
                                      : gelex::GeneticMode::D)}
                            : std::string{"{add|dom}"};

    auto output_pattern
        = config.loco
              ? fmt::format(
                    "{}.{}.chr{{01..{:02d}}}.{{bin|id}}",
                    config.out,
                    task_pattern,
                    ranges.size())
              : fmt::format("{}.{}.{{bin|id}}", config.out, task_pattern);

    cli::GrmReporter::show_files_written(
        ranges.size() * modes.size() * 2,
        std::filesystem::absolute(std::filesystem::path(config.out))
            .parent_path()
            .string(),
        output_pattern);

    return 0;
}
