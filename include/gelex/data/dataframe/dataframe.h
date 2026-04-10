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
#include <initializer_list>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <fmt/format.h>

#include "gelex/data/dataframe/column.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/dataframe/key_type.h"
#include "gelex/exception.h"
#include "gelex/infra/string_hash.h"

namespace gelex::df
{

namespace detail
{
template <KeyType Key>
class DataFrameReader;
}  // namespace detail

template <KeyType Key>
class DataFrame
{
   public:
    DataFrame(const DataFrame&) = delete;
    auto operator=(const DataFrame&) -> DataFrame& = delete;
    DataFrame(DataFrame&&) = default;
    auto operator=(DataFrame&&) -> DataFrame& = default;
    ~DataFrame() = default;

    // to get column by name;
    template <typename Self>
    auto&& operator[](this Self&& self, std::string_view name);
    // get column by index;
    template <typename Self>
    auto&& col(this Self&& self, std::size_t index);

    auto names() const -> std::vector<std::string>
    {
        return columns_
               | std::views::transform([](const auto& c)
                                       { return std::string(c.name()); })
               | std::ranges::to<std::vector>();
    }
    auto cols() const -> std::size_t { return columns_.size(); }
    auto rows() const -> std::size_t { return index_.size(); }
    template <typename Self>
    auto&& index(this Self&& self)
    {
        return std::forward<Self>(self).index_;
    }

    // check if column exists
    auto contains(std::string_view name) const -> bool
    {
        return col_lookup_.contains(name);
    }

    // rename columns, length of new_names must be equal to columns
    auto rename(std::span<const std::string> new_names) -> void;

    auto clone() const -> DataFrame;

    // gather rows by gelex::df::Index
    auto gather(const Index<Key>& target) -> void;

    template <ValueType T = double>
        requires std::is_arithmetic_v<T>
    auto to_mat() const -> Eigen::MatrixX<T>;

    template <ValueType T = double>
        requires std::is_arithmetic_v<T>
    auto to_mat(std::span<const std::size_t> col_indices) const
        -> Eigen::MatrixX<T>;

    template <ValueType T = double>
        requires std::is_arithmetic_v<T>
    auto to_mat(std::span<const std::string_view> col_names) const
        -> Eigen::MatrixX<T>;

   private:
    DataFrame() = default;
    friend detail::DataFrameReader<Key>;

    // internal gather directly by index
    auto gather(std::span<const std::size_t> indices) -> void;

    auto push_back(Column col) -> void;
    auto update_lookup(std::size_t index, std::string_view new_name) -> void;
    auto set_name(std::size_t index, std::string_view new_name) -> void;

    Index<Key> index_;
    std::vector<Column> columns_;
    std::unordered_map<
        std::string,
        std::size_t,
        infra::TransparentHash<std::string>,
        infra::TransparentEqual<std::string>>
        col_lookup_;
};

template <KeyType Key>
auto intersect_inplace(std::span<DataFrame<Key>* const> dfs) -> void;

template <KeyType Key>
auto intersect_inplace(std::initializer_list<DataFrame<Key>* const> dfs) -> void
{
    intersect_inplace(std::span{dfs.begin(), dfs.size()});
}

// --- Implementation ---

template <KeyType Key>
template <typename Self>
auto&& DataFrame<Key>::operator[](this Self&& self, std::string_view name)
{
    auto it = self.col_lookup_.find(name);
    if (it == self.col_lookup_.end())
    {
        throw GelexException(fmt::format("column not found: '{}'", name));
    }
    return std::forward<Self>(self).col(it->second);
}

template <KeyType Key>
template <typename Self>
auto&& DataFrame<Key>::col(this Self&& self, std::size_t index)
{
    if (index >= self.columns_.size())
    {
        throw GelexException(
            fmt::format("column index out of range: {}", index));
    }
    return std::forward<Self>(self).columns_[index];
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
    decltype(col_lookup_) new_lookup;
    for (std::size_t i = 0; i < new_names.size(); ++i)
    {
        auto [it, ok] = new_lookup.emplace(new_names[i], i);
        if (!ok)
        {
            throw GelexException(
                fmt::format("duplicate column name: '{}'", new_names[i]));
        }
    }
    col_lookup_ = std::move(new_lookup);
    for (auto&& [name, col] : std::views::zip(new_names, columns_))
    {
        col.rename(name);
    }
}

template <KeyType Key>
auto DataFrame<Key>::clone() const -> DataFrame
{
    DataFrame result;
    result.index_ = index_;
    result.columns_ = columns_;
    result.col_lookup_ = col_lookup_;
    return result;
}

template <KeyType Key>
auto DataFrame<Key>::gather(const Index<Key>& target) -> void
{
    auto pos = target.keys()
               | std::views::transform([this](const auto& k)
                                       { return index_.at(k); })
               | std::ranges::to<std::vector>();
    for (auto& c : columns_)
    {
        c.gather(pos);
    }
    index_ = target;
}

template <KeyType Key>
template <ValueType T>
    requires std::is_arithmetic_v<T>
auto DataFrame<Key>::to_mat() const -> Eigen::MatrixX<T>
{
    auto c = static_cast<Eigen::Index>(columns_.size());
    Eigen::MatrixX<T> mat(static_cast<Eigen::Index>(rows()), c);
    for (Eigen::Index j = 0; j < c; ++j)
    {
        mat.col(j) = col(static_cast<std::size_t>(j)).template to_map<T>();
    }
    return mat;
}

template <KeyType Key>
template <ValueType T>
    requires std::is_arithmetic_v<T>
auto DataFrame<Key>::to_mat(std::span<const std::size_t> col_indices) const
    -> Eigen::MatrixX<T>
{
    auto c = static_cast<Eigen::Index>(col_indices.size());
    Eigen::MatrixX<T> mat(static_cast<Eigen::Index>(rows()), c);
    for (Eigen::Index j = 0; j < c; ++j)
    {
        mat.col(j) = col(col_indices[static_cast<std::size_t>(j)])
                         .template to_map<T>();
    }
    return mat;
}

template <KeyType Key>
template <ValueType T>
    requires std::is_arithmetic_v<T>
auto DataFrame<Key>::to_mat(std::span<const std::string_view> col_names) const
    -> Eigen::MatrixX<T>
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
auto DataFrame<Key>::update_lookup(std::size_t index, std::string_view new_name)
    -> void
{
    auto [it, inserted] = col_lookup_.emplace(new_name, index);
    if (!inserted)
    {
        throw GelexException(
            fmt::format("duplicate column name: '{}'", new_name));
    }
}

template <KeyType Key>
auto DataFrame<Key>::push_back(Column col) -> void
{
    auto idx = columns_.size();
    update_lookup(idx, col.name());
    columns_.push_back(std::move(col));
}

template <KeyType Key>
auto DataFrame<Key>::set_name(std::size_t index, std::string_view new_name)
    -> void
{
    update_lookup(index, new_name);
    columns_[index].rename(new_name);
}

template <KeyType Key>
auto intersect_inplace(std::span<DataFrame<Key>* const> dfs) -> void
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
    auto common = intersect(std::span<const Index<Key>* const>{idx_ptrs});
    for (auto* df : dfs)
    {
        df->gather(common);
    }
}

}  // namespace gelex::df

#endif  // GELEX_DATA_DATAFRAME_DATAFRAME_H
