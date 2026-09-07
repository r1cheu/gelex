// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_DATA_SAMPLE_ID_H_
#define GELEX_DATA_SAMPLE_ID_H_

#include <string_view>
#include <utility>

#include "gelex/exception.h"

namespace gelex
{

inline constexpr char sample_id_separator = '\x1F';

inline auto split_sample_id(std::string_view sample_id)
    -> std::pair<std::string_view, std::string_view>
{
    if (sample_id.empty())
    {
        throw GelexException("sample ID cannot be empty");
    }

    auto separator_pos = sample_id.find(sample_id_separator);
    if (separator_pos == std::string_view::npos)
    {
        throw GelexException(
            "sample ID is not in the canonical FID<US>IID format");
    }

    auto fid = sample_id.substr(0, separator_pos);
    auto iid = sample_id.substr(separator_pos + 1);
    if (fid.empty() || iid.empty())
    {
        throw GelexException("sample ID must contain non-empty FID and IID");
    }

    return {fid, iid};
}

}  // namespace gelex

#endif  // GELEX_DATA_SAMPLE_ID_H_
