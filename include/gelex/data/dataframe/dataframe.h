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

#ifndef GELEX_DATA_DATAFRAME_DATAFRAME_H
#define GELEX_DATA_DATAFRAME_DATAFRAME_H

#include <cstddef>
#include <format>
#include <initializer_list>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "gelex/data/dataframe/column.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/dataframe/key_type.h"
#include "gelex/data/dataframe/string_hash.h"
#include "gelex/exception.h"

namespace gelex::df
{

template <KeyType Key>
class DataFrameReader;

template <KeyType Key>
class DataFrame
{
   public:
    ~DataFrame() = default;
    DataFrame(const DataFrame&) = delete;
    auto operator=(const DataFrame&) -> DataFrame& = delete;
    DataFrame(DataFrame&&) = default;
    auto operator=(DataFrame&&) -> DataFrame& = default;

    template <typename Self>
    auto&& col(this Self&& self, std::string_view name)
    {
        auto it = self.col_lookup_.find(name);
        if (it == self.col_lookup_.end())
        {
            throw InvalidInputException(
                std::format("column not found: '{}'", name));
        }
        return std::forward<Self>(self).col(it->second);
    }
    template <typename Self>
    auto&& col(this Self&& self, std::size_t index)
    {
        if (index >= self.columns_.size())
        {
            throw InvalidInputException(
                std::format("column index out of range: {}", index));
        }
        return std::forward<Self>(self).columns_[index];
    }

    auto names() const -> std::span<const std::string> { return names_; }
    auto cols() const -> std::size_t { return columns_.size(); }
    auto contains(std::string_view name) const -> bool
    {
        return col_lookup_.contains(name);
    }

    template <typename Self>
    auto&& index(this Self&& self)
    {
        return std::forward<Self>(self).index_;
    }
    auto rows() const -> std::size_t { return index_.size(); }
    auto row_position(const Key& key) const -> std::size_t
    {
        return index_[key];
    }

    auto clone() const -> DataFrame;
    auto gather(std::span<const std::size_t> indices) -> void;

    static auto intersect(std::span<DataFrame* const> dfs) -> void;
    static auto intersect(std::initializer_list<DataFrame* const> dfs) -> void;

    template <ValueType T = double>
        requires std::is_arithmetic_v<T>
    auto to_mat(std::span<const std::size_t> col_indices = {}) const
        -> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    template <ValueType T = double>
        requires std::is_arithmetic_v<T>
    auto to_mat(std::span<const std::string_view> col_names) const
        -> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

   private:
    DataFrame() = default;
    friend DataFrameReader<Key>;

    Index<Key> index_;
    std::vector<std::string> names_;
    std::vector<Column> columns_;
    std::unordered_map<
        std::string,
        std::size_t,
        TransparentHash<std::string>,
        TransparentEqual<std::string>>
        col_lookup_;
};

// --- Implementation ---

template <KeyType Key>
auto DataFrame<Key>::clone() const -> DataFrame
{
    DataFrame result;
    result.index_ = index_;
    result.names_ = names_;
    result.columns_ = columns_;
    result.col_lookup_ = col_lookup_;
    return result;
}

template <KeyType Key>
auto DataFrame<Key>::gather(std::span<const std::size_t> indices) -> void
{
    for (auto& col : columns_)
    {
        col.gather(indices);
    }
    index_.gather(indices);
}

template <KeyType Key>
auto DataFrame<Key>::intersect(std::span<DataFrame* const> dfs) -> void
{
    if (dfs.empty())
    {
        return;
    }

    std::vector<const Index<Key>*> idx_ptrs;
    idx_ptrs.reserve(dfs.size());
    for (const auto* df : dfs)
    {
        idx_ptrs.push_back(&df->index_);
    }
    auto positions = Index<Key>::intersect(idx_ptrs);

    for (auto&& [df, pos] : std::views::zip(dfs, positions))
    {
        df->gather(pos);
    }
}

template <KeyType Key>
auto DataFrame<Key>::intersect(std::initializer_list<DataFrame* const> dfs)
    -> void
{
    intersect(std::span{dfs.begin(), dfs.size()});
}

template <KeyType Key>
template <ValueType T>
    requires std::is_arithmetic_v<T>
auto DataFrame<Key>::to_mat(std::span<const std::size_t> col_indices) const
    -> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
{
    auto c = col_indices.empty()
                 ? static_cast<Eigen::Index>(columns_.size())
                 : static_cast<Eigen::Index>(col_indices.size());
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat(
        static_cast<Eigen::Index>(rows()), c);
    for (Eigen::Index j = 0; j < c; ++j)
    {
        auto idx = col_indices.empty()
                       ? static_cast<std::size_t>(j)
                       : col_indices[static_cast<std::size_t>(j)];
        auto data = col(idx).template as<T>();
        mat.col(j) = Eigen::Map<const Eigen::Vector<T, Eigen::Dynamic>>(
            data.data(), static_cast<Eigen::Index>(data.size()));
    }
    return mat;
}

template <KeyType Key>
template <ValueType T>
    requires std::is_arithmetic_v<T>
auto DataFrame<Key>::to_mat(std::span<const std::string_view> col_names) const
    -> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
{
    std::vector<std::size_t> indices;
    indices.reserve(col_names.size());
    for (auto name : col_names)
    {
        auto it = col_lookup_.find(name);
        if (it == col_lookup_.end())
        {
            throw InvalidInputException(
                std::format("column not found: '{}'", name));
        }
        indices.push_back(it->second);
    }
    return to_mat<T>(std::span<const std::size_t>{indices});
}

}  // namespace gelex::df

#endif  // GELEX_DATA_DATAFRAME_DATAFRAME_H
