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

#ifndef GELEX_DATA_DATAFRAME_COLUMN_H
#define GELEX_DATA_DATAFRAME_COLUMN_H

#include <fmt/format.h>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/dataframe/key_type.h"
#include "gelex/exception.h"

namespace gelex::df
{

enum class ColumnType : std::uint8_t
{
    Int,
    Float,
    Double,
    String
};

using ColumnStorage = std::variant<
    std::monostate,
    std::vector<int32_t>,
    std::vector<float>,
    std::vector<double>,
    std::vector<std::string>>;

class Column
{
   public:
    explicit Column(std::string_view name) : name_(name) {}

    template <ValueType T>
    Column(std::string_view name, std::vector<T> data)
        : name_(name), storage_(std::move(data))
    {
    }

    template <ValueType T>
    auto push_back(T&& value) -> void;

    auto name() const -> std::string_view { return name_; }
    auto rename(std::string_view new_name) -> void { name_ = new_name; }

    auto size() const -> std::size_t;

    // reorders rows to match the given index order; indices may be in any order
    auto gather(std::span<const std::size_t> indices) -> void;

    template <ValueType T>
    auto as() -> std::span<T>;

    template <ValueType T>
    auto as() const -> std::span<const T>;

    template <ValueType T>
    auto take() && -> std::vector<T>;

    template <ValueType T>
        requires std::is_arithmetic_v<T>
    auto to_map() const -> Eigen::Map<const Eigen::Vector<T, Eigen::Dynamic>>;

    template <ValueType T>
        requires std::is_arithmetic_v<T>
    auto to_mat() const -> Eigen::Vector<T, Eigen::Dynamic>;

   private:
    template <ValueType T, typename Storage>
    static auto checked(Storage& storage, std::string_view col_name) -> auto&
    {
        auto* vec = std::get_if<std::vector<T>>(&storage);
        if (!vec)
        {
            throw GelexException(
                fmt::format("column '{}': type mismatch", col_name));
        }
        return *vec;
    }

    std::string name_;
    ColumnStorage storage_;
};

// --- Implementation ---

template <ValueType T>
auto Column::push_back(T&& value) -> void
{
    using Raw = std::remove_cvref_t<T>;
    if (std::holds_alternative<std::monostate>(storage_))
    {
        storage_ = std::vector<Raw>{std::forward<T>(value)};
        return;
    }
    auto* vec = std::get_if<std::vector<Raw>>(&storage_);
    if (!vec)
    {
        throw GelexException(
            fmt::format("column '{}': push_back type mismatch", name_));
    }
    vec->push_back(std::forward<T>(value));
}

template <ValueType T>
auto Column::as() -> std::span<T>
{
    return checked<T>(storage_, name_);
}

template <ValueType T>
auto Column::as() const -> std::span<const T>
{
    return checked<T>(storage_, name_);
}

template <ValueType T>
auto Column::take() && -> std::vector<T>
{
    return std::move(checked<T>(storage_, name_));
}

template <ValueType T>
    requires std::is_arithmetic_v<T>
auto Column::to_map() const
    -> Eigen::Map<const Eigen::Vector<T, Eigen::Dynamic>>
{
    const auto& vec = checked<T>(storage_, name_);
    return {vec.data(), static_cast<Eigen::Index>(vec.size())};
}

template <ValueType T>
    requires std::is_arithmetic_v<T>
auto Column::to_mat() const -> Eigen::Vector<T, Eigen::Dynamic>
{
    return to_map<T>();
}

}  // namespace gelex::df

#endif  // GELEX_DATA_DATAFRAME_COLUMN_H
