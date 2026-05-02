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

#include "gelex/data/dataframe/column.h"

#include <cassert>
#include <cstddef>
#include <span>
#include <type_traits>
#include <variant>

namespace gelex::dataframe
{

auto Column::size() const -> std::size_t
{
    return std::visit(
        [](const auto& v) -> std::size_t
        {
            if constexpr (std::is_same_v<
                              std::decay_t<decltype(v)>,
                              std::monostate>)
            {
                return 0;
            }
            else
            {
                return v.size();
            }
        },
        storage_);
}

auto Column::gather(std::span<const std::size_t> indices) -> void
{
    std::visit(
        [&indices](auto& v)
        {
            if constexpr (!std::is_same_v<
                              std::decay_t<decltype(v)>,
                              std::monostate>)
            {
                std::decay_t<decltype(v)> tmp;
                tmp.reserve(indices.size());
                for (auto i : indices)
                {
                    assert(i < v.size());
                    tmp.push_back(v[i]);
                }
                v = std::move(tmp);
            }
        },
        storage_);
}

}  // namespace gelex::dataframe
