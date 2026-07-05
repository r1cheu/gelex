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
#include <optional>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/bayes/labels.h"
#include "gelex/data/genotype.h"
#include "gelex/data/genotype_reader.h"
#include "gelex/exception.h"
#include "gelex/io/binary_reader.h"
#include "gelex/post/diagnostic.h"
#include "gelex/post/fixed.h"
#include "gelex/post/genetic.h"
#include "gelex/post/genetic_variance.h"
#include "gelex/post/genetic_variance_kernel.h"
#include "gelex/post/heritability.h"
#include "gelex/post/random.h"
#include "gelex/post/residual.h"
#include "gelex/types/genetic_mode.h"
#include "reporter.h"

namespace
{

auto parse_genetic_mode(std::string_view value) -> gelex::GeneticMode
{
    for (const auto& [mode, name] : gelex::GENETIC_MODE_NAMES)
    {
        if (value == name)
        {
            return mode;
        }
    }
    throw gelex::GelexException(
        fmt::format("unknown genetic mode in samples: {}", value));
}

auto check_consistency(const std::vector<gelex::io::BinaryReader>& readers)
    -> bool
{
    if (readers.size() <= 1)
    {
        return !readers.empty();
    }

    const auto reference = readers.front().section_paths();
    return std::ranges::all_of(
        readers | std::views::drop(1),
        [&](const gelex::io::BinaryReader& reader)
        { return reader.section_paths() == reference; });
}

auto process_gebv_variance(
    std::vector<gelex::io::BinaryReader>& readers,
    const std::optional<std::string>& gfile,
    double hdpi_threshold) -> std::vector<gelex::ParameterDiag>
{
    if (!gfile)
    {
        return {};
    }

    std::vector<gelex::ParameterDiag> diags;

    std::vector<gelex::Genotype> genotype_storages;
    std::vector<gelex::GeneticInput> genetic_inputs;
    const auto& ref = readers.front();
    if (!ref.contains("genetic/modes"))
    {
        return {};
    }

    const auto sample_modes = ref.to_strings("genetic/modes");
    genotype_storages.reserve(sample_modes.size());
    for (const auto [index, mode_name] : std::views::enumerate(sample_modes))
    {
        const auto kind = parse_genetic_mode(mode_name);
        auto coeff_path = fmt::format("genetic/{}/coeffs", index);
        if (!readers.front().contains(coeff_path))
        {
            continue;
        }
        auto gbin_path = fmt::format(
            "{}.{}.gbin", *gfile, gelex::bayes::to_file_suffix(kind));

        genotype_storages.emplace_back(
            gelex::genotype::GenotypeReader::read(gbin_path, kind));
        genetic_inputs.push_back({&genotype_storages.back(), kind});
    }

    if (!genetic_inputs.empty())
    {
        auto [gebv_diags, genetic_variances]
            = gelex::
                  GeneticVariancePosteriorProcessor{readers, genetic_inputs, hdpi_threshold}
                      .process();
        diags.append_range(std::move(gebv_diags));

        std::vector<gelex::GeneticMode> active_kinds;
        active_kinds.reserve(genetic_inputs.size());
        for (const auto& input : genetic_inputs)
        {
            active_kinds.push_back(input.kind);
        }
        diags.append_range(
            gelex::HeritabilityPosteriorProcessor{
                readers, hdpi_threshold, genetic_variances, active_kinds}
                .process());
    }

    return diags;
}

}  // namespace

auto post_execute(const cli::PostConfig& config) -> int
{
    cli::PostReporter reporter;
    const auto& in_prefixes = config.in;

    std::vector<gelex::io::BinaryReader> readers;
    readers.reserve(in_prefixes.size());
    for (const auto& prefix : in_prefixes)
    {
        readers.emplace_back(fmt::format("{}.samples", prefix));
    }

    if (!check_consistency(readers))
    {
        throw gelex::GelexException(
            "Inconsistent section paths across sample files");
    }

    const auto hdpi_threshold = config.hdpi;
    std::span<const gelex::io::BinaryReader> reader_span{readers};

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

    auto gebv_diags
        = process_gebv_variance(readers, config.gfile, hdpi_threshold);

    std::vector<gelex::ParameterDiag> diags;
    diags.append_range(std::move(fixed_diags));
    diags.append_range(std::move(random_diags));
    diags.append_range(std::move(gebv_diags));
    diags.append_range(std::move(genetic_diags));
    diags.append_range(std::move(residual_diags));

    reporter.show_diagnostics(diags, hdpi_threshold);

    return 0;
}
