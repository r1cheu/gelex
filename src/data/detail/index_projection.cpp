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

#include "gelex/data/detail/index_projection.h"

#include <cstddef>
#include <fmt/format.h>

#include "gelex/exception.h"

namespace gelex::detail
{
IndexProjection::IndexProjection(
    const DataFrameIndex<std::string>& source_index,
    const DataFrameIndex<std::string>& target_index)
    : source_size_{static_cast<index_type>(source_index.size())},
      source_to_target_(static_cast<std::size_t>(source_size_), npos),
      is_identity_{source_index.size() == target_index.size()}
{
    target_to_source_.reserve(target_index.size());

    for (std::size_t target_pos = 0; target_pos < target_index.size();
         ++target_pos)
    {
        const auto& sample_id = target_index.keys()[target_pos];

        const auto source_pos
            = static_cast<index_type>(source_index.at(sample_id));

        target_to_source_.push_back(source_pos);

        if (source_to_target_[static_cast<std::size_t>(source_pos)] != npos)
        {
            throw GelexException{
                fmt::format("duplicated projected sample: {}", sample_id)};
        }

        source_to_target_[static_cast<std::size_t>(source_pos)]
            = static_cast<index_type>(target_pos);

        if (static_cast<std::size_t>(source_pos) != target_pos)
        {
            is_identity_ = false;
        }
    }
}

}  // namespace gelex::detail
