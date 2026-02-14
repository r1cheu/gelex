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

#ifndef GELEX_DATA_DATAFRAME_POLICY_H_
#define GELEX_DATA_DATAFRAME_POLICY_H_

#include <cstdint>
#include <format>
#include <functional>
#include <string>
#include <string_view>
#include <unordered_set>

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

enum class MissingValueAction : uint8_t
{
    Throw,
    UseTypeDefault,
    SkipRow,
};

struct TransparentStringHash
{
    using is_transparent = void;

    [[nodiscard]] auto operator()(const std::string& value) const noexcept
        -> std::size_t
    {
        return std::hash<std::string>{}(value);
    }

    [[nodiscard]] auto operator()(std::string_view value) const noexcept
        -> std::size_t
    {
        return std::hash<std::string_view>{}(value);
    }
};

using MissingTokenSet
    = std::unordered_set<std::string, TransparentStringHash, std::equal_to<>>;

inline const MissingTokenSet kDefaultMissingTokens
    = {"", "NA", "NaN", "nan", "null", "NULL", "."};

struct DataFrameLoadPolicy
{
    MissingTokenSet missing_tokens = kDefaultMissingTokens;
    MissingValueAction missing_value_action = MissingValueAction::SkipRow;
};

}  // namespace gelex

#endif  // GELEX_DATA_DATAFRAME_POLICY_H_
