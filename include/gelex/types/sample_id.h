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

#ifndef GELEX_TYPES_SAMPLE_ID_H_
#define GELEX_TYPES_SAMPLE_ID_H_

#include <format>
#include <string>
#include <string_view>
#include <utility>

#include "gelex/exception.h"

namespace gelex
{

inline constexpr char kSampleIdSeparator = '\x1F';

inline auto make_sample_id(std::string_view fid, std::string_view iid)
    -> std::string
{
    if (fid.empty())
    {
        throw ArgumentValidationException("FID cannot be empty");
    }
    if (iid.empty())
    {
        throw ArgumentValidationException("IID cannot be empty");
    }

    return std::format("{}{}{}", fid, kSampleIdSeparator, iid);
}

inline auto split_sample_id(std::string_view sample_id)
    -> std::pair<std::string_view, std::string_view>
{
    if (sample_id.empty())
    {
        throw ArgumentValidationException("sample ID cannot be empty");
    }

    auto separator_pos = sample_id.find(kSampleIdSeparator);
    if (separator_pos == std::string_view::npos)
    {
        throw ArgumentValidationException(
            "sample ID is not in the canonical FID<US>IID format");
    }

    auto fid = sample_id.substr(0, separator_pos);
    auto iid = sample_id.substr(separator_pos + 1);
    if (fid.empty() || iid.empty())
    {
        throw ArgumentValidationException(
            "sample ID must contain non-empty FID and IID");
    }

    return {fid, iid};
}

}  // namespace gelex

#endif  // GELEX_TYPES_SAMPLE_ID_H_
