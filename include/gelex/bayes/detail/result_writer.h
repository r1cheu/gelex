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

#ifndef GELEX_BAYES_DETAIL_RESULT_WRITER_H_
#define GELEX_BAYES_DETAIL_RESULT_WRITER_H_

#include <Eigen/Core>
#include <cstddef>
#include <fmt/format.h>
#include <ranges>
#include <string>
#include <string_view>

#include "gelex/bayes/basic_result.h"
#include "gelex/bayes/genetic/result.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/detail/text_writer.h"

namespace gelex::detail
{

inline auto format_result_value(double value) -> std::string
{
    return fmt::format("{:.8e}", value);
}

inline auto write_parameter_rows(
    TextWriter& writer,
    const CoefficientPosteriorResult& result) -> void
{
    const auto& statistics = result.statistics();
    for (const auto [index, column_name] :
         result.column_names() | std::views::enumerate)
    {
        const auto statistics_index = static_cast<Eigen::Index>(index);
        writer.write(
            fmt::format(
                "{}\t{}\t{}\t{}\t{}",
                result.identifier(),
                index,
                column_name,
                format_result_value(statistics.mean(statistics_index)),
                format_result_value(statistics.stddev(statistics_index))));
    }
}

inline auto write_summary_row(
    TextWriter& writer,
    std::string_view identifier,
    std::size_t index,
    double mean,
    double stddev) -> void
{
    writer.write(
        fmt::format(
            "{}\t{}\t{}\t{}",
            identifier,
            index,
            format_result_value(mean),
            format_result_value(stddev)));
}

// Marker-sized posteriors stay in the binary trace: the summary file holds
// one row per scalar term, not one per marker.
inline auto write_summary_rows(
    TextWriter& /*writer*/,
    const EmptyPosteriorResult& /*result*/) -> void
{
}

inline auto write_summary_rows(
    TextWriter& writer,
    const ScalarPosteriorResult& result) -> void
{
    const auto& statistics = result.statistics();
    write_summary_row(
        writer, result.identifier(), 0, statistics.mean, statistics.stddev);
}

inline auto write_summary_rows(
    TextWriter& writer,
    const VectorPosteriorResult& result) -> void
{
    const auto& statistics = result.statistics();
    for (const auto [index, mean] : statistics.mean | std::views::enumerate)
    {
        const auto statistics_index = static_cast<Eigen::Index>(index);
        write_summary_row(
            writer,
            result.identifier(),
            static_cast<std::size_t>(index),
            mean,
            statistics.stddev(statistics_index));
    }
}

template <VarianceLayout Kind>
auto write_family_summary_rows(
    TextWriter& writer,
    const GaussianPosteriorResult<Kind>& result) -> void
{
    write_summary_rows(writer, result.variance);
}

inline auto write_family_summary_rows(
    TextWriter& writer,
    const HalfNormalPosteriorResult& result) -> void
{
    write_summary_rows(writer, result.variance);
    write_summary_rows(writer, result.probit_coefficients);
}

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
auto write_family_summary_rows(
    TextWriter& writer,
    const SpikeSlabPosteriorResult<Kind, WeightUpdate>& result) -> void
{
    write_summary_rows(writer, result.variance);
    write_summary_rows(writer, result.probability);
}

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
auto write_family_summary_rows(
    TextWriter& writer,
    const ScaledMixturePosteriorResult<ClassCount, WeightUpdate>& result)
    -> void
{
    write_summary_rows(writer, result.variance);
    write_summary_rows(writer, result.probabilities);
    write_summary_rows(writer, result.component_explained_variance);
}

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
auto write_family_summary_rows(
    TextWriter& writer,
    const JointSpikeSlabPosteriorResult<ClassCount, WeightUpdate>& result)
    -> void
{
    write_summary_rows(writer, result.probabilities);
    write_summary_rows(writer, result.component_explained_variance);
}

template <GeneticModeSet Modes, typename... Results>
auto write_genetic_summary_rows(
    TextWriter& writer,
    const ModeValues<Modes, Results...>& result) -> void
{
    result.for_each([&]<GeneticMode Mode>(const auto& mode_result)
                    { write_family_summary_rows(writer, mode_result); });
}

template <typename ModeValuesType, typename JointResult>
auto write_genetic_summary_rows(
    TextWriter& writer,
    const JointModeValues<ModeValuesType, JointResult>& result) -> void
{
    write_genetic_summary_rows(writer, result.mode_values());
    write_family_summary_rows(writer, result.joint());
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_DETAIL_RESULT_WRITER_H_
