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

#ifndef GELEX_DATA_DATAFRAME_H_
#define GELEX_DATA_DATAFRAME_H_

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <format>
#include <span>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "gelex/data/column.h"
#include "gelex/data/dataframe_policy.h"
#include "gelex/exception.h"
#include "gelex/internal/dataframe/dataframe_loader.h"
#include "gelex/internal/dataframe/index_intersector.h"

namespace gelex::detail
{

template <typename T>
class DataFrameLoader;

template <typename T>
struct IndexIntersector;

}  // namespace gelex::detail

namespace gelex
{

template <typename T>
class DataFrame
{
   public:
    DataFrame();

    static auto read(
        const std::filesystem::path& path,
        const DataFrameLoadPolicy& policy = {}) -> DataFrame;

    auto intersect_index_inplace(std::span<const std::string> keys) -> void;

    [[nodiscard]] auto nrows() const -> size_t;

    [[nodiscard]] auto ncols() const -> size_t;

    [[nodiscard]] auto index_column() const -> const Column<std::string>&;

    [[nodiscard]] auto column(size_t i) const -> const Column<T>&;

   private:
    auto initialize_columns(std::vector<std::string> names) -> void;

    auto rebuild_index_map() -> void;

    Column<std::string> index_;
    std::vector<Column<T>> columns_;
    std::unordered_map<std::string, size_t> index_map_;

    friend class detail::DataFrameLoader<T>;
    friend struct detail::IndexIntersector<T>;
};

template <typename T>
auto compute_common_index_keys(std::span<const DataFrame<T>* const> frames)
    -> std::vector<std::string>
{
    if (frames.empty())
    {
        return {};
    }

    if (frames[0] == nullptr)
    {
        throw ArgumentValidationException("frames[0] cannot be null");
    }

    std::vector<std::string> common = frames[0]->index_column().data();

    for (size_t i = 1; i < frames.size(); ++i)
    {
        const DataFrame<T>* frame = frames[i];
        if (frame == nullptr)
        {
            throw ArgumentValidationException(
                std::format("frames[{}] cannot be null", i));
        }

        std::unordered_set<std::string> next_keys;
        const auto& next_index = frame->index_column().data();
        next_keys.reserve(next_index.size());
        next_keys.insert(next_index.begin(), next_index.end());

        common.erase(
            std::remove_if(
                common.begin(),
                common.end(),
                [&next_keys](const std::string& key)
                { return !next_keys.contains(key); }),
            common.end());

        if (common.empty())
        {
            break;
        }
    }

    return common;
}

}  // namespace gelex

namespace gelex
{

template <typename T>
DataFrame<T>::DataFrame() : index_("index")
{
}

template <typename T>
auto DataFrame<T>::read(
    const std::filesystem::path& path,
    const DataFrameLoadPolicy& policy) -> DataFrame
{
    return detail::DataFrameLoader<T>::load(path, policy);
}

template <typename T>
auto DataFrame<T>::intersect_index_inplace(std::span<const std::string> keys)
    -> void
{
    detail::IndexIntersector<T>::apply(*this, keys);
}

template <typename T>
auto DataFrame<T>::nrows() const -> size_t
{
    return index_.size();
}

template <typename T>
auto DataFrame<T>::ncols() const -> size_t
{
    return columns_.size() + 1;
}

template <typename T>
auto DataFrame<T>::index_column() const -> const Column<std::string>&
{
    return index_;
}

template <typename T>
auto DataFrame<T>::column(size_t i) const -> const Column<T>&
{
    if (i >= columns_.size())
    {
        throw ColumnRangeException(std::format("column {} is out of range", i));
    }

    return columns_[i];
}

template <typename T>
auto DataFrame<T>::initialize_columns(std::vector<std::string> names) -> void
{
    columns_.clear();
    columns_.reserve(names.size());
    for (auto& name : names)
    {
        columns_.emplace_back(std::move(name));
    }
}

template <typename T>
auto DataFrame<T>::rebuild_index_map() -> void
{
    index_map_.clear();
    index_map_.reserve(index_.size());

    const auto& idx = index_.data();
    for (size_t i = 0; i < idx.size(); ++i)
    {
        auto [_, inserted] = index_map_.emplace(idx[i], i);
        if (!inserted)
        {
            throw InvalidOperationException("duplicated index key");
        }
    }
}

}  // namespace gelex

#endif  // GELEX_DATA_DATAFRAME_H_
