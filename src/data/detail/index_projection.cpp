// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

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
