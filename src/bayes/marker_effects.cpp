// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/bayes/marker_effects.h"

#include <Eigen/Core>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <iterator>
#include <string>
#include <string_view>

#include "gelex/bayes/genetic/marker_effect_traits.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/exception.h"
#include "gelex/io/detail/atomic_output_stream.h"

namespace gelex
{

auto write_marker_effects(
    std::string_view path,
    const bayes::GeneticDesign& design,
    const MarkerEffectTable& table) -> void
{
    if (table.rows() != design.cols())
    {
        throw GelexException(
            fmt::format(
                "write_marker_effects: table has {} rows but the design has "
                "{} markers",
                table.rows(),
                design.cols()));
    }
    const auto& metadata = design.marker_metadata();
    const auto snp = metadata.index().keys();
    const auto chromosome = metadata["CHR"].as<std::string>();
    const auto position = metadata["BP"].as<std::int32_t>();
    const auto a1 = metadata["A1"].as<std::string>();
    const auto a2 = metadata["A2"].as<std::string>();
    const auto& a1_frequency = design.a1_frequency();

    detail::AtomicOutputStream file{std::filesystem::path{path}};
    std::string line{"CHR\tSNP\tBP\tA1\tA2\tA1FREQ"};
    for (const auto& name : table.names())
    {
        fmt::format_to(std::back_inserter(line), "\t{}", name);
    }
    line.push_back('\n');
    file.write(line);

    for (Eigen::Index marker = 0; marker < table.rows(); ++marker)
    {
        const auto row = static_cast<std::size_t>(marker);
        line.clear();
        fmt::format_to(
            std::back_inserter(line),
            "{}\t{}\t{}\t{}\t{}\t{:.10g}",
            chromosome[row],
            snp[row],
            position[row],
            a1[row],
            a2[row],
            a1_frequency(marker));
        for (const auto& column : table.columns())
        {
            fmt::format_to(
                std::back_inserter(line), "\t{:.10g}", column(marker));
        }
        line.push_back('\n');
        file.write(line);
    }
    file.commit();
}

}  // namespace gelex
