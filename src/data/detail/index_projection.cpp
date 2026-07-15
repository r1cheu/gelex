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
#include <vector>

#include "gelex/exception.h"

namespace gelex::detail
{
IndexProjection::IndexProjection(
    const DataFrameIndex<std::string>& source_index,
    const DataFrameIndex<std::string>& target_index)
    : source_size_{static_cast<index_type>(source_index.size())}
{
    target_to_source_.reserve(target_index.size());

    std::vector<bool> projected(static_cast<std::size_t>(source_size_), false);

    for (const auto& sample_id : target_index.keys())
    {
        const auto source_pos
            = static_cast<index_type>(source_index.at(sample_id));

        if (projected[static_cast<std::size_t>(source_pos)])
        {
            throw GelexException{
                fmt::format("duplicated projected sample: {}", sample_id)};
        }
        projected[static_cast<std::size_t>(source_pos)] = true;

        target_to_source_.push_back(source_pos);
    }
}

}  // namespace gelex::detail
