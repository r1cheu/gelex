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

#ifndef GELEX_DATA_COLUMN_H_
#define GELEX_DATA_COLUMN_H_

#include <cstddef>
#include <format>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "gelex/exception.h"

namespace gelex
{

template <typename T>
class Column
{
   public:
    explicit Column(std::string name = "") : name_(std::move(name)) {}

    auto append(T value) -> void { data_.push_back(std::move(value)); }

    auto gather_inplace(std::span<const size_t> keep_rows) -> void
    {
        std::vector<T> filtered;
        filtered.reserve(keep_rows.size());

        for (size_t row : keep_rows)
        {
            if (row >= data_.size())
            {
                throw ColumnRangeException(
                    std::format("row {} is out of range", row));
            }
            filtered.emplace_back(std::move(data_[row]));
        }

        data_ = std::move(filtered);
    }

    [[nodiscard]] auto size() const -> size_t { return data_.size(); }

    [[nodiscard]] auto name() const -> const std::string& { return name_; }

    [[nodiscard]] auto data() const -> const std::vector<T>& { return data_; }

    [[nodiscard]] auto data() -> std::vector<T>& { return data_; }

   private:
    std::string name_;
    std::vector<T> data_;
};

}  // namespace gelex

#endif  // GELEX_DATA_COLUMN_H_
