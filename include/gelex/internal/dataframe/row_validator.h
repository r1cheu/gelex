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

#ifndef GELEX_INTERNAL_DATAFRAME_ROW_VALIDATOR_H_
#define GELEX_INTERNAL_DATAFRAME_ROW_VALIDATOR_H_

#include <cstddef>
#include <format>
#include <string_view>

#include "gelex/exception.h"

namespace gelex::detail
{

inline auto validate_expected_columns(
    size_t actual,
    size_t expected,
    size_t line_number) -> void
{
    if (actual != expected)
    {
        throw InconsistentColumnCountException(
            std::format(
                "line {} has {} columns, expected {}",
                line_number,
                actual,
                expected));
    }
}

inline auto validate_minimum_columns(
    size_t actual,
    size_t minimum,
    size_t line_number) -> void
{
    if (actual < minimum)
    {
        throw InconsistentColumnCountException(
            std::format(
                "line {} has {} columns, expected {}",
                line_number,
                actual,
                minimum));
    }
}

inline auto validate_non_empty_token(
    std::string_view token,
    size_t line_number,
    std::string_view token_name) -> void
{
    if (token.empty())
    {
        throw DataParseException(
            std::format("line {} has missing {}", line_number, token_name));
    }
}

}  // namespace gelex::detail

#endif  // GELEX_INTERNAL_DATAFRAME_ROW_VALIDATOR_H_
