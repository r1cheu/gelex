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

#ifndef GELEX_BAYES_GENETIC_DETAIL_SUMMARY_SUPPORT_H_
#define GELEX_BAYES_GENETIC_DETAIL_SUMMARY_SUPPORT_H_

#include <Eigen/Core>
#include <cstddef>
#include <fmt/format.h>
#include <ranges>
#include <string>
#include <string_view>

#include "gelex/bayes/basic_result.h"
#include "gelex/io/detail/text_writer.h"

namespace gelex::detail
{

inline auto format_result_value(double value) -> std::string
{
    return fmt::format("{:.8e}", value);
}

inline auto write_parameter_rows(
    TextWriter& writer,
    const CoefficientResult& result) -> void
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
    const EmptyResult& /*result*/) -> void
{
}

inline auto write_summary_rows(TextWriter& writer, const ScalarResult& result)
    -> void
{
    const auto& statistics = result.statistics();
    write_summary_row(
        writer, result.identifier(), 0, statistics.mean, statistics.stddev);
}

inline auto write_summary_rows(TextWriter& writer, const VectorResult& result)
    -> void
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

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_DETAIL_SUMMARY_SUPPORT_H_
