// Copyright 2026 RuLei Chen
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "gelex/bayes/marker_covariate.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cstdint>
#include <fmt/format.h>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "gelex/data/dataframe/constants.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/exception.h"

namespace
{
template <gelex::ValueType T>
auto validate_marker_metadata_column(
    const gelex::DataFrame<std::string>& frame,
    const gelex::DataFrame<std::string>& marker_metadata,
    std::string_view column_name) -> void
{
    const auto found_values = frame[column_name].as<T>();
    const auto expected_values = marker_metadata[column_name].as<T>();
    const auto [found, expected]
        = std::ranges::mismatch(found_values, expected_values);
    if (found == found_values.end())
    {
        return;
    }

    const auto marker = *std::ranges::next(
        marker_metadata.index().keys().begin(),
        std::ranges::distance(found_values.begin(), found));
    throw gelex::GelexException(
        fmt::format(
            "marker '{}' {} mismatch: expected '{}', found '{}'",
            marker,
            column_name,
            *expected,
            *found));
}
}  // namespace

namespace gelex::bayes
{
auto make_marker_covariate(
    DataFrame<std::string> frame,
    const DataFrame<std::string>& marker_metadata) -> MarkerCovariate
{
    constexpr std::array<std::string_view, 4> metadata_names{
        "CHR", "BP", "A1", "A2"};
    const auto frame_names = frame.names();
    if (frame_names.size() <= metadata_names.size())
    {
        throw GelexException(
            "marker covariate requires at least one annotation column");
    }
    if (!std::ranges::equal(
            frame_names | std::views::take(metadata_names.size()),
            metadata_names))
    {
        throw GelexException(
            "marker covariate metadata columns must be CHR, BP, A1, A2");
    }

    frame.gather(marker_metadata.index());
    validate_marker_metadata_column<std::string>(frame, marker_metadata, "CHR");
    validate_marker_metadata_column<std::int32_t>(frame, marker_metadata, "BP");
    validate_marker_metadata_column<std::string>(frame, marker_metadata, "A1");
    validate_marker_metadata_column<std::string>(frame, marker_metadata, "A2");

    const auto annotation_count = frame_names.size() - metadata_names.size();
    Eigen::MatrixXd values(
        static_cast<Eigen::Index>(annotation_count + 1),
        static_cast<Eigen::Index>(frame.rows()));
    values.row(0).setOnes();

    std::vector<std::string> annotation_names;
    annotation_names.reserve(annotation_count + 1);
    annotation_names.emplace_back(intercept_name);
    for (const auto [index, column_name] : std::views::enumerate(
             frame_names | std::views::drop(metadata_names.size())))
    {
        annotation_names.push_back(column_name);
        values.row(static_cast<Eigen::Index>(index + 1))
            = frame[column_name].to_map<double>().transpose();
    }

    if (!values.allFinite())
    {
        throw GelexException("marker covariate values must be finite");
    }

    return MarkerCovariate(std::move(annotation_names), std::move(values));
}
}  // namespace gelex::bayes
