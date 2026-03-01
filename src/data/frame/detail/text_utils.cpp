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

#include "gelex/data/frame/detail/text_utils.h"

#include <cstddef>
#include <string_view>
#include <vector>

namespace gelex::detail
{

auto split_line_preserve_empty(std::string_view line, char delimiter)
    -> std::vector<std::string_view>
{
    std::vector<std::string_view> tokens;
    size_t start = 0;

    while (true)
    {
        size_t pos = line.find(delimiter, start);
        if (pos == std::string_view::npos)
        {
            tokens.emplace_back(line.substr(start));
            break;
        }
        tokens.emplace_back(line.substr(start, pos - start));
        start = pos + 1;
    }

    return tokens;
}

auto detect_delimiter(std::string_view probe) -> char
{
    if (probe.find('\t') != std::string_view::npos)
    {
        return '\t';
    }
    if (probe.find(',') != std::string_view::npos)
    {
        return ',';
    }
    return ' ';
}

}  // namespace gelex::detail
