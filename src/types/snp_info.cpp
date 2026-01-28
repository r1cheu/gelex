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

#include "gelex/types/snp_info.h"

namespace gelex
{

SnpEffects::SnpEffects(size_t initial_capacity)
{
    if (initial_capacity > 0)
    {
        snp_meta_.reserve(initial_capacity);
        snp_index_map_.reserve(initial_capacity);

        additive_data_.reserve(initial_capacity);
        frequencies_data_.reserve(initial_capacity);
    }
}

auto SnpEffects::emplace_meta(SnpMeta meta) -> void
{
    snp_index_map_.emplace(
        meta.id, static_cast<Eigen::Index>(snp_meta_.size()));
    snp_meta_.emplace_back(std::move(meta));
}

auto SnpEffects::find_index(std::string_view snp_id) const
    -> std::optional<Eigen::Index>
{
    auto it = snp_index_map_.find(std::string(snp_id));
    if (it != snp_index_map_.end())
    {
        return it->second;
    }
    return std::nullopt;
}

auto SnpEffects::operator[](std::string_view snp_id) -> SnpMeta*
{
    auto idx = find_index(snp_id);
    if (idx)
    {
        return &snp_meta_[*idx];
    }
    return nullptr;
}

auto SnpEffects::operator[](std::string_view snp_id) const -> const SnpMeta*
{
    auto idx = find_index(snp_id);
    if (idx)
    {
        return &snp_meta_[*idx];
    }
    return nullptr;
}

auto SnpEffects::shrink_to_fit() -> void
{
    snp_meta_.shrink_to_fit();
    // unordered_map 没有 shrink_to_fit，通常忽略

    additive_data_.shrink_to_fit();
    frequencies_data_.shrink_to_fit();
    dominance_data_.shrink_to_fit();
}

auto SnpEffects::clear() -> void
{
    snp_meta_.clear();
    snp_index_map_.clear();
    additive_data_.clear();
    frequencies_data_.clear();
    dominance_data_.clear();
}

}  // namespace gelex
