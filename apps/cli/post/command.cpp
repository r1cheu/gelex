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

#include <algorithm>
#include <fmt/format.h>
#include <ranges>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "gelex/exception.h"
#include "gelex/io/binary_reader.h"
#include "gelex/post/diagnostic.h"
#include "gelex/post/fixed.h"
#include "gelex/post/genetic.h"
#include "gelex/post/random.h"
#include "gelex/post/residual.h"

#include "reporter.h"

namespace
{

auto check_consistency(const std::vector<gelex::BinaryReader>& readers) -> bool
{
    if (readers.size() <= 1)
    {
        return !readers.empty();
    }

    const auto reference = readers.front().section_paths();
    return std::ranges::all_of(
        readers | std::views::drop(1),
        [&](const gelex::BinaryReader& reader)
        { return reader.section_paths() == reference; });
}

}  // namespace

auto post_execute(const cli::PostConfig& config) -> int
{
    cli::PostReporter reporter;
    const auto& in_prefixes = config.in;

    std::vector<gelex::BinaryReader> readers;
    readers.reserve(in_prefixes.size());
    for (const auto& prefix : in_prefixes)
    {
        readers.emplace_back(fmt::format("{}.draws", prefix));
    }

    if (!check_consistency(readers))
    {
        throw gelex::GelexException(
            "Inconsistent section paths across sample files");
    }

    const auto hdpi_threshold = config.hdpi;
    std::span<const gelex::BinaryReader> reader_span{readers};

    auto fixed_diags
        = gelex::FixedPosteriorProcessor{reader_span, hdpi_threshold}.process();
    auto random_diags
        = gelex::RandomPosteriorProcessor{reader_span, hdpi_threshold}
              .process();
    auto genetic_diags
        = gelex::GeneticPosteriorProcessor{reader_span, hdpi_threshold}
              .process();
    auto residual_diags
        = gelex::ResidualPosteriorProcessor{reader_span, hdpi_threshold}
              .process();

    std::vector<gelex::ParameterDiag> diags;
    diags.append_range(std::move(fixed_diags));
    diags.append_range(std::move(random_diags));
    diags.append_range(std::move(genetic_diags));
    diags.append_range(std::move(residual_diags));

    reporter.show_diagnostics(diags, hdpi_threshold);

    return 0;
}
