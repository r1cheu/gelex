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

#include "gelex/data/sample_manager.h"

#include <algorithm>
#include <filesystem>
#include <ranges>
#include <utility>
#include <vector>

#include "loader/fam_loader.h"

namespace gelex
{

SampleManager::SampleManager(
    const std::filesystem::path& fam_path,
    bool iid_only)
{
    detail::FamLoader fam_loader(fam_path, iid_only);

    common_ids_ = std::move(fam_loader).take_ids();

    std::ranges::sort(common_ids_);

    const auto [first, last] = std::ranges::unique(common_ids_);
    common_ids_.erase(first, last);
}

void SampleManager::intersect(std::span<const std::string> ids)
{
    if (common_ids_.empty())
    {
        return;
    }
    if (ids.empty())
    {
        common_ids_.clear();
        return;
    }

    std::vector<std::string_view> sorted_input(ids.begin(), ids.end());
    std::ranges::sort(sorted_input);

    std::erase_if(
        common_ids_,
        [&sorted_input](const std::string& s)
        {
            return !std::ranges::binary_search(sorted_input, s, std::less<>{});
        });
}

void SampleManager::finalize()
{
    common_id_map_.clear();
    common_id_map_.reserve(common_ids_.size());

    for (auto&& [idx, id] : common_ids_ | std::views::enumerate)
    {
        common_id_map_.emplace(id, static_cast<Eigen::Index>(idx));
    }
}

}  // namespace gelex
