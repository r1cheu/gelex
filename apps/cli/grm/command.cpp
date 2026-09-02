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

#include <Eigen/Core>
#include <cstddef>
#include <fmt/format.h>
#include <string>
#include <vector>

#include "gelex/data/bed.h"
#include "gelex/data/grm/grm.h"
#include "gelex/data/grm/io.h"
#include "gelex/data/marker_range.h"
#include "gelex/genetic_mode.h"

#include "cli/formatter.h"
#include "cli/report_printer.h"
#include "cli/runtime.h"
#include "reporter.h"

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
        ranges = gelex::chromosome_ranges(bed.bim());
    }
    else
    {
        ranges.push_back({std::string{}, 0, bed.num_snps()});
    }

    cli::printer().block(cli::section("GRM Computation:"));
    gelex::GrmBuilder builder(bed, modes, method, chunk_size, observer);
    builder.build(
        ranges,
        [&](const gelex::GrmMatrix& matrix)
        {
            const std::string name = matrix.label.empty()
                                         ? fmt::format("{}", matrix.mode)
                                         : fmt::format(
                                               "{}.chr{:02d}",
                                               matrix.mode,
                                               std::stoi(matrix.label));
            gelex::write_grm(
                fmt::format("{}.{}", config.out, name), matrix.grm, sample_ids);
        });

    reporter.finish_progress();

    const std::string task_pattern
        = modes.size() == 1 ? fmt::format(
                                  "{}",
                                  modes.contains(gelex::GeneticMode::A)
                                      ? gelex::GeneticMode::A
                                      : gelex::GeneticMode::D)
                            : std::string{"{A,D}"};

    auto suffix_pattern = config.loco
                              ? fmt::format(
                                    ".{}.chr{{01..{:02d}}}.{{bin,id}}",
                                    task_pattern,
                                    ranges.size())
                              : fmt::format(".{}.{{bin,id}}", task_pattern);

    cli::printer().block(
        cli::results_saved(
            config.out, fmt::format("{}, .log", suffix_pattern)));

    return 0;
}
