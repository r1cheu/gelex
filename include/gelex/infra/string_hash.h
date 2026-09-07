// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_INFRA_STRING_HASH_H_
#define GELEX_INFRA_STRING_HASH_H_

#include <cstddef>
#include <functional>
#include <string>
#include <string_view>
#include <unordered_set>

namespace gelex
{

template <typename Key>
struct TransparentHash
{
    using is_transparent = void;
    auto operator()(const Key& k) const -> std::size_t
    {
        return std::hash<Key>{}(k);
    }
};

template <>
struct TransparentHash<std::string>
{
    using is_transparent = void;
    auto operator()(std::string_view sv) const -> std::size_t
    {
        return std::hash<std::string_view>{}(sv);
    }
};

template <typename Key>
struct TransparentEqual
{
    using is_transparent = void;
    auto operator()(const Key& a, const Key& b) const -> bool { return a == b; }
};

template <>
struct TransparentEqual<std::string>
{
    using is_transparent = void;
    auto operator()(std::string_view a, std::string_view b) const -> bool
    {
        return a == b;
    }
};

using StringSet = std::unordered_set<
    std::string,
    TransparentHash<std::string>,
    TransparentEqual<std::string>>;

}  // namespace gelex

#endif  // GELEX_INFRA_STRING_HASH_H_
