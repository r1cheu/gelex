// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/data/dataframe/column.h"

#include <cassert>
#include <cstddef>
#include <span>
#include <type_traits>
#include <variant>

namespace gelex
{

auto Column::size() const -> std::size_t
{
    return std::visit(
        [](const auto& v) -> std::size_t
        {
            if constexpr (
                std::is_same_v<std::decay_t<decltype(v)>, std::monostate>)
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

}  // namespace gelex
