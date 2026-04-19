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

#ifndef GELEX_DATA_DATAFRAME_INDEX_H
#define GELEX_DATA_DATAFRAME_INDEX_H

#include <fmt/format.h>
#include <algorithm>
#include <concepts>
#include <cstddef>
#include <initializer_list>
#include <ranges>
#include <span>
#include <unordered_map>
#include <vector>

#include "gelex/data/dataframe/key_type.h"
#include "gelex/exception.h"
#include "gelex/infra/string_hash.h"

namespace gelex::df
{

template <KeyType Key>
class Index
{
   public:
    Index() = default;
    explicit Index(std::vector<Key> keys);

    template <typename K>
        requires std::convertible_to<K, Key>
    auto push_back(K&& key) -> void;

    auto at(const Key& key) const -> std::size_t;
    auto contains(const Key& key) const -> bool
    {
        return lookup_.contains(key);
    }
    auto size() const -> std::size_t { return keys_.size(); }

    auto keys() const -> std::span<const Key> { return keys_; }
    auto take_keys() && -> std::vector<Key> { return std::move(keys_); }

    auto gather(const Index<Key>& target) -> void;
    auto gather(std::span<const std::size_t> indices) -> void;

   private:
    auto rebuild_lookup() -> void;

    template <KeyType K>
    friend auto intersect(std::span<const Index<K>* const> indices) -> Index<K>;

    std::vector<Key> keys_;
    std::unordered_map<
        Key,
        std::size_t,
        infra::TransparentHash<Key>,
        infra::TransparentEqual<Key>>
        lookup_;
};

template <KeyType Key>
auto intersect(std::span<const Index<Key>* const> indices) -> Index<Key>;

template <KeyType Key>
auto intersect(std::initializer_list<const Index<Key>* const> indices)
    -> Index<Key>
{
    return intersect(std::span{indices.begin(), indices.size()});
}

// --- Implementation ---

template <KeyType Key>
Index<Key>::Index(std::vector<Key> keys) : keys_(std::move(keys))
{
    for (std::size_t i = 0; i < keys_.size(); ++i)
    {
        const auto& key = keys_[i];
        if (lookup_.contains(key))
        {
            throw GelexException(fmt::format("duplicate index: {}", key));
        }
        lookup_[key] = i;
    }
}

template <KeyType Key>
template <typename K>
    requires std::convertible_to<K, Key>
auto Index<Key>::push_back(K&& key) -> void
{
    if (lookup_.contains(key))
    {
        throw GelexException(fmt::format("duplicate index: {}", key));
    }
    keys_.push_back(std::forward<K>(key));
    lookup_[keys_.back()] = keys_.size() - 1;
}

template <KeyType Key>
auto Index<Key>::at(const Key& key) const -> std::size_t
{
    auto it = lookup_.find(key);
    if (it == lookup_.end())
    {
        throw GelexException(fmt::format("index not found: {}", key));
    }
    return it->second;
}

template <KeyType Key>
auto Index<Key>::gather(const Index<Key>& target) -> void
{
    auto pos = target.keys()
               | std::views::transform([this](const auto& k) { return at(k); })
               | std::ranges::to<std::vector>();
    gather(pos);
}

template <KeyType Key>
auto Index<Key>::gather(std::span<const std::size_t> indices) -> void
{
    std::vector<Key> buf;
    buf.reserve(indices.size());
    for (const auto i : indices)
    {
        buf.push_back(keys_[i]);
    }
    keys_ = std::move(buf);
    rebuild_lookup();
}

template <KeyType Key>
auto Index<Key>::rebuild_lookup() -> void
{
    lookup_.clear();
    lookup_.reserve(keys_.size());
    for (std::size_t i = 0; i < keys_.size(); ++i)
    {
        lookup_[keys_[i]] = i;
    }
}

template <KeyType Key>
auto intersect(std::span<const Index<Key>* const> indices) -> Index<Key>
{
    if (indices.empty())
    {
        return Index<Key>{};
    }

    const auto* smallest = *std::ranges::min_element(
        indices,
        [](const auto* a, const auto* b) { return a->size() < b->size(); });

    std::vector<Key> common;
    common.reserve(smallest->size());
    for (const auto& key : smallest->keys_)
    {
        bool in_all = std::ranges::all_of(
            indices,
            [&key, smallest](const auto* idx)
            { return idx == smallest || idx->lookup_.contains(key); });
        if (in_all)
        {
            common.push_back(key);
        }
    }
    std::ranges::sort(common);
    return Index<Key>(std::move(common));
}

}  // namespace gelex::df

#endif  // GELEX_DATA_DATAFRAME_INDEX_H
