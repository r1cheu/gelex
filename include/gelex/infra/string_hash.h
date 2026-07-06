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
