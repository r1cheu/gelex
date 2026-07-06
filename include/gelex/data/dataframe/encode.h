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

#ifndef GELEX_DATA_DATAFRAME_ENCODE_H_
#define GELEX_DATA_DATAFRAME_ENCODE_H_

#include <fmt/format.h>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/dataframe/column.h"
#include "gelex/exception.h"

namespace gelex
{

template <typename Scalar = double>
    requires std::is_floating_point_v<Scalar>
struct EncodedResult
{
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> data;
    std::string name;
    std::vector<std::string> level_names;
};

struct LevelMismatch
{
    std::vector<std::string> missing_in_data;
    std::vector<std::string> missing_in_levels;

    [[nodiscard]] auto ok() const -> bool
    {
        return missing_in_data.empty() && missing_in_levels.empty();
    }
};

namespace detail
{

[[nodiscard]] auto collect_unique_sorted(std::span<const std::string> values)
    -> std::vector<std::string>;

[[nodiscard]] auto make_level_names(
    std::string_view col_name,
    std::span<const std::string> levels) -> std::vector<std::string>;

[[nodiscard]] auto has_duplicates(std::span<const std::string> levels) -> bool;

template <typename Scalar>
[[nodiscard]] auto encode_matrix(
    std::span<const std::string> values,
    std::span<const std::string> levels)
    -> Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
{
    const auto n_rows = static_cast<Eigen::Index>(values.size());
    const auto n_cols = static_cast<Eigen::Index>(levels.size());

    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> mat(n_rows, n_cols);
    mat.setZero();

    std::unordered_map<std::string, Eigen::Index> level_idx;
    level_idx.reserve(levels.size());
    for (Eigen::Index i = 0; i < n_cols; ++i)
    {
        level_idx.emplace(levels[static_cast<std::size_t>(i)], i);
    }

    for (Eigen::Index row = 0; row < n_rows; ++row)
    {
        const auto it = level_idx.find(values[static_cast<std::size_t>(row)]);
        if (it != level_idx.end())
        {
            mat(row, it->second) = Scalar{1};
        }
    }

    return mat;
}

}  // namespace detail

[[nodiscard]] inline auto collect_levels(const Column& col)
    -> std::vector<std::string>
{
    return detail::collect_unique_sorted(col.as<std::string>());
}

[[nodiscard]] auto check_levels(
    const Column& col,
    std::span<const std::string> levels) -> LevelMismatch;

// --- Implementation ---

template <typename Scalar = double>
    requires std::is_floating_point_v<Scalar>
[[nodiscard]] auto encode(
    const Column& col,
    std::span<const std::string> levels) -> EncodedResult<Scalar>
{
    if (levels.empty())
    {
        throw GelexException(
            fmt::format(
                "encode: levels must not be empty for column '{}'",
                col.name()));
    }
    if (detail::has_duplicates(levels))
    {
        throw GelexException(
            fmt::format(
                "encode: levels contain duplicates for column '{}'",
                col.name()));
    }

    auto values = col.as<std::string>();
    return EncodedResult<Scalar>{
        .data = detail::encode_matrix<Scalar>(values, levels),
        .name = std::string(col.name()),
        .level_names = detail::make_level_names(col.name(), levels)};
}

template <typename Scalar = double>
    requires std::is_floating_point_v<Scalar>
[[nodiscard]] auto one_hot_encode(const Column& col) -> EncodedResult<Scalar>
{
    auto levels = collect_levels(col);
    if (levels.size() < 2)
    {
        throw GelexException(
            fmt::format(
                "encoding requires at least 2 levels, column '{}' has {}",
                col.name(),
                levels.size()));
    }
    return encode<Scalar>(col, levels);
}

template <typename Scalar = double>
    requires std::is_floating_point_v<Scalar>
[[nodiscard]] auto dummy_encode(const Column& col) -> EncodedResult<Scalar>
{
    auto levels = collect_levels(col);
    if (levels.size() < 2)
    {
        throw GelexException(
            fmt::format(
                "encoding requires at least 2 levels, column '{}' has {}",
                col.name(),
                levels.size()));
    }
    return encode<Scalar>(col, std::span<const std::string>(levels).subspan(1));
}

}  // namespace gelex

#endif  // GELEX_DATA_DATAFRAME_ENCODE_H_
