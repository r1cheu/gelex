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

#ifndef GELEX_DATA_FRAME_DETAIL_INDEX_INTERSECTOR_H_
#define GELEX_DATA_FRAME_DETAIL_INDEX_INTERSECTOR_H_

#include <span>
#include <string>
#include <unordered_set>
#include <vector>

#include "gelex/exception.h"

namespace gelex
{

template <typename T>
class DataFrame;

namespace detail
{

template <typename T>
struct IndexIntersector
{
    static auto apply(DataFrame<T>& frame, std::span<const std::string> keys)
        -> void
    {
        std::unordered_set<std::string> key_set;
        key_set.reserve(keys.size());

        for (const auto& key : keys)
        {
            auto [_, inserted] = key_set.insert(key);
            if (!inserted)
            {
                throw InvalidOperationException(
                    "duplicate key in intersect input");
            }
        }

        std::vector<size_t> keep_rows;
        keep_rows.reserve(keys.size());

        for (const auto& key : keys)
        {
            auto it = frame.index_map_.find(key);
            if (it != frame.index_map_.end())
            {
                keep_rows.push_back(it->second);
            }
        }

        frame.index_.gather_inplace(keep_rows);
        for (auto& column : frame.columns_)
        {
            column.gather_inplace(keep_rows);
        }
        frame.rebuild_index_map();
    }
};

}  // namespace detail

}  // namespace gelex

#endif  // GELEX_DATA_FRAME_DETAIL_INDEX_INTERSECTOR_H_
