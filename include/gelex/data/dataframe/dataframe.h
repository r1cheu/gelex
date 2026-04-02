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

#include <fmt/format.h>
#include <cstddef>
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
#include "gelex/exception.h"
#include "gelex/infra/string_hash.h"

namespace gelex::df
{

template <KeyType Key>
class DataFrameReader;

template <KeyType Key>
class DataFrame
{
   public:
    DataFrame(const DataFrame&) = delete;
    auto operator=(const DataFrame&) -> DataFrame& = delete;
    DataFrame(DataFrame&&) = default;
    auto operator=(DataFrame&&) -> DataFrame& = default;
    ~DataFrame() = default;

    template <typename Self>
    auto&& operator[](this Self&& self, std::string_view name)
    {
        auto it = self.col_lookup_.find(name);
        if (it == self.col_lookup_.end())
        {
            throw GelexException(fmt::format("column not found: '{}'", name));
        }
        return std::forward<Self>(self).col(it->second);
    }
    template <typename Self>
    auto&& col(this Self&& self, std::size_t index)
    {
        if (index >= self.columns_.size())
        {
            throw GelexException(
                fmt::format("column index out of range: {}", index));
        }
        return std::forward<Self>(self).columns_[index];
    }

    auto names() const -> std::span<const std::string> { return names_; }
    auto cols() const -> std::size_t { return columns_.size(); }
    auto contains(std::string_view name) const -> bool
    {
        return col_lookup_.contains(name);
    }

    auto rename(std::span<const std::string> new_names) -> void;

    template <typename Self>
    auto&& index(this Self&& self)
    {
        return std::forward<Self>(self).index_;
    }
    auto rows() const -> std::size_t { return index_.size(); }

    auto clone() const -> DataFrame;
    auto gather(std::span<const std::size_t> indices) -> void;

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

    auto set_name(std::size_t index, std::string_view new_name) -> void;

    Index<Key> index_;
    std::vector<std::string> names_;
    std::vector<Column> columns_;
    std::unordered_map<
        std::string,
        std::size_t,
        infra::TransparentHash<std::string>,
        infra::TransparentEqual<std::string>>
        col_lookup_;
};

template <KeyType Key>
auto intersect(std::span<DataFrame<Key>* const> dfs) -> void;

template <KeyType Key>
auto intersect(std::initializer_list<DataFrame<Key>* const> dfs) -> void
{
    intersect(std::span{dfs.begin(), dfs.size()});
}

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
auto DataFrame<Key>::set_name(std::size_t index, std::string_view new_name)
    -> void
{
    std::string name(new_name);
    auto [it, inserted] = col_lookup_.emplace(name, index);
    if (!inserted)
    {
        throw GelexException(
            fmt::format("duplicate column name: '{}'", new_name));
    }
    names_[index] = it->first;
    columns_[index].rename(it->first);
}

template <KeyType Key>
auto DataFrame<Key>::rename(std::span<const std::string> new_names) -> void
{
    if (new_names.size() != columns_.size())
    {
        throw GelexException(
            fmt::format(
                "rename: {} names given but DataFrame has {} columns",
                new_names.size(),
                columns_.size()));
    }
    col_lookup_.clear();
    for (std::size_t i = 0; i < new_names.size(); ++i)
    {
        set_name(i, new_names[i]);
    }
}

template <KeyType Key>
auto DataFrame<Key>::gather(std::span<const std::size_t> indices) -> void
{
    for (auto& c : columns_)
    {
        c.gather(indices);
    }
    index_.gather(indices);
}

template <KeyType Key>
auto intersect(std::span<DataFrame<Key>* const> dfs) -> void
{
    if (dfs.empty())
    {
        return;
    }

    std::vector<const Index<Key>*> idx_ptrs;
    idx_ptrs.reserve(dfs.size());
    for (const auto* df : dfs)
    {
        idx_ptrs.push_back(&df->index());
    }
    auto positions = intersect(std::span<const Index<Key>* const>{idx_ptrs});

    for (auto&& [df, pos] : std::views::zip(dfs, positions))
    {
        df->gather(pos);
    }
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
            throw GelexException(fmt::format("column not found: '{}'", name));
        }
        indices.push_back(it->second);
    }
    return to_mat<T>(std::span<const std::size_t>{indices});
}

}  // namespace gelex::df

#endif  // GELEX_DATA_DATAFRAME_DATAFRAME_H
