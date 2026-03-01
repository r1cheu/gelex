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

#ifndef GELEX_DATA_FRAME_DETAIL_VALUE_PARSER_H_
#define GELEX_DATA_FRAME_DETAIL_VALUE_PARSER_H_

#include <charconv>
#include <concepts>
#include <format>
#include <string>
#include <string_view>

#include "gelex/exception.h"

namespace gelex::detail
{

template <typename T>
struct ValueParser
{
    static auto parse(std::string_view token) -> T
    {
        static_assert(
            std::integral<T> || std::floating_point<T>,
            "ValueParser<T> only supports arithmetic types by default");

        T value{};
        const char* first = token.data();
        const char* last = token.data() + token.size();
        auto [ptr, ec] = std::from_chars(first, last, value);
        if (ec == std::errc() && ptr == last)
        {
            return value;
        }

        throw NumberParseException(
            std::format("failed to parse '{}' as number", token));
    }
};

template <>
struct ValueParser<std::string>
{
    static auto parse(std::string_view token) -> std::string
    {
        return std::string(token);
    }
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_FRAME_DETAIL_VALUE_PARSER_H_
