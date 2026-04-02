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

#include "gelex/data/dataframe/encode.h"

#include <fmt/format.h>
#include <algorithm>
#include <span>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

#include "gelex/data/dataframe/column.h"
#include "gelex/data/dataframe/constants.h"

namespace gelex::df
{

namespace detail
{

auto collect_unique_sorted(std::span<const std::string> values)
    -> std::vector<std::string>
{
    std::unordered_set<std::string> seen(values.begin(), values.end());
    std::vector<std::string> levels(seen.begin(), seen.end());
    std::sort(levels.begin(), levels.end());
    return levels;
}

auto make_level_names(
    std::string_view col_name,
    std::span<const std::string> levels) -> std::vector<std::string>
{
    std::vector<std::string> names;
    names.reserve(levels.size());
    for (const auto& level : levels)
    {
        names.push_back(fmt::format("{}{}{}", col_name, kSeparator, level));
    }
    return names;
}

auto has_duplicates(std::span<const std::string> levels) -> bool
{
    std::unordered_set<std::string> seen;
    seen.reserve(levels.size());
    for (const auto& lv : levels)
    {
        if (!seen.insert(lv).second)
        {
            return true;
        }
    }
    return false;
}

}  // namespace detail

auto check_levels(const Column& col, std::span<const std::string> levels)
    -> LevelMismatch
{
    const auto data_levels
        = detail::collect_unique_sorted(col.as<std::string>());

    const std::unordered_set<std::string> data_set(
        data_levels.begin(), data_levels.end());
    const std::unordered_set<std::string> level_set(
        levels.begin(), levels.end());

    LevelMismatch mismatch;

    for (const auto& lv : levels)
    {
        if (!data_set.contains(lv))
        {
            mismatch.missing_in_data.push_back(lv);
        }
    }
    for (const auto& lv : data_levels)
    {
        if (!level_set.contains(lv))
        {
            mismatch.missing_in_levels.push_back(lv);
        }
    }

    std::sort(mismatch.missing_in_data.begin(), mismatch.missing_in_data.end());
    std::sort(
        mismatch.missing_in_levels.begin(), mismatch.missing_in_levels.end());

    return mismatch;
}

}  // namespace gelex::df
