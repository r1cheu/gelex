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

#ifndef GELEX_BAYES_RESULT_IO_H_
#define GELEX_BAYES_RESULT_IO_H_

#include <Eigen/Core>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <iterator>
#include <string>
#include <string_view>
#include <type_traits>

#include "gelex/bayes/detail/result_writer.h"
#include "gelex/bayes/genetic/result.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/result.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/detail/text_writer.h"

namespace gelex
{

template <typename GeneticPrior>
auto write_params(
    const BayesResult<GeneticPrior>& result,
    std::string_view prefix) -> void
{
    detail::TextWriter writer{fmt::format("{}.params", prefix)};
    writer.write_header({"term", "index", "column_name", "mean", "stddev"});

    detail::write_parameter_rows(writer, result.fixed());
    for (const auto& random : result.random())
    {
        detail::write_parameter_rows(writer, random.coefficients());
    }
}

template <typename GeneticPrior>
auto write_summary(
    const BayesResult<GeneticPrior>& result,
    std::string_view prefix) -> void
{
    detail::TextWriter writer{fmt::format("{}.summary", prefix)};
    writer.write_header({"term", "index", "mean", "stddev"});

    for (const auto& random : result.random())
    {
        detail::write_summary_rows(writer, random.variance());
        detail::write_summary_rows(writer, random.explained_variance());
    }

    detail::write_genetic_summary_rows(writer, result.genetic_parameters());
    detail::write_summary_rows(writer, result.residual());

    const auto& variance_summary = result.variance_summary();
    const auto write_mode_summary = [&writer]<GeneticMode>(const auto& summary)
    { detail::write_summary_rows(writer, summary); };
    variance_summary.explained_variances().for_each(write_mode_summary);
    variance_summary.heritabilities().for_each(write_mode_summary);
    detail::write_summary_rows(
        writer, variance_summary.total_explained_variance());
    detail::write_summary_rows(writer, variance_summary.total_heritability());
}

template <typename GeneticPrior>
auto write_snpeff(
    const BayesResult<GeneticPrior>& result,
    const bayes::GeneticDesign& genetic_design,
    std::string_view prefix) -> void
{
    const auto& marker_effects = result.marker_effects();

    std::string header{"CHR\tSNP\tBP\tA1\tA2\tA1FREQ"};
    marker_effects.for_each(
        [&header]<GeneticMode Mode>(const auto& mode_effects)
        {
            fmt::format_to(
                std::back_inserter(header),
                "\tBETA_{}\tSE_{}\tPVE_{}",
                Mode,
                Mode,
                Mode);
            using PipResult = std::remove_cvref_t<decltype(mode_effects.pip())>;
            if constexpr (std::same_as<PipResult, MarkerPipResult>)
            {
                fmt::format_to(std::back_inserter(header), "\tPIP_{}", Mode);
            }
        });
    if constexpr (requires { marker_effects.joint(); })
    {
        header += "\tPVE";
        using PipResult
            = std::remove_cvref_t<decltype(marker_effects.joint().pip())>;
        if constexpr (std::same_as<PipResult, MarkerPipResult>)
        {
            header += "\tPIP";
        }
    }

    detail::TextWriter writer{fmt::format("{}.snpeff", prefix)};
    writer.write(header);

    const auto& marker_metadata = genetic_design.marker_metadata();
    const auto marker_ids = marker_metadata.index().keys();
    const auto chromosomes = marker_metadata["chrom"].as<std::string>();
    const auto positions = marker_metadata["pos"].as<std::int32_t>();
    const auto a1 = marker_metadata["A1"].as<std::string>();
    const auto a2 = marker_metadata["A2"].as<std::string>();
    const auto& a1_frequency = genetic_design.a1_frequency();

    std::string line;
    line.reserve(128 + (marker_effects.modes.size() * 64));
    for (Eigen::Index marker = 0; marker < genetic_design.cols(); ++marker)
    {
        const auto row = static_cast<std::size_t>(marker);
        line.clear();
        fmt::format_to(
            std::back_inserter(line),
            "{}\t{}\t{}\t{}\t{}\t{:.8e}",
            chromosomes[row],
            marker_ids[row],
            positions[row],
            a1[row],
            a2[row],
            a1_frequency(marker));

        marker_effects.for_each(
            [&line, marker]<GeneticMode>(const auto& mode_effects)
            {
                const auto& statistics
                    = mode_effects.coefficients().statistics();
                fmt::format_to(
                    std::back_inserter(line),
                    "\t{:.8e}\t{:.8e}\t{:.8e}",
                    statistics.mean(marker),
                    statistics.stddev(marker),
                    mode_effects.pve().values()(marker));
                using PipResult
                    = std::remove_cvref_t<decltype(mode_effects.pip())>;
                if constexpr (std::same_as<PipResult, MarkerPipResult>)
                {
                    fmt::format_to(
                        std::back_inserter(line),
                        "\t{:.8e}",
                        mode_effects.pip().probabilities()(marker));
                }
            });
        if constexpr (requires { marker_effects.joint(); })
        {
            const auto& joint = marker_effects.joint();
            fmt::format_to(
                std::back_inserter(line),
                "\t{:.8e}",
                joint.pve().values()(marker));
            using PipResult = std::remove_cvref_t<decltype(joint.pip())>;
            if constexpr (std::same_as<PipResult, MarkerPipResult>)
            {
                fmt::format_to(
                    std::back_inserter(line),
                    "\t{:.8e}",
                    joint.pip().probabilities()(marker));
            }
        }
        writer.write(line);
    }
}

}  // namespace gelex

#endif  // GELEX_BAYES_RESULT_IO_H_
