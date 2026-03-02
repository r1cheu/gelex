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

#include "gelex/io/parser.h"

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <ranges>
#include <vector>

#include "gelex/types/sample_id.h"

namespace gelex::detail
{

size_t count_total_lines(const std::filesystem::path& path)
{
    constexpr size_t buffer_size = 128 * 1024;

    std::vector<char> buffer(buffer_size);

    auto file = open_file<std::ifstream>(path, std::ios::in | std::ios::binary);

    file.rdbuf()->pubsetbuf(buffer.data(), buffer_size);

    size_t line_count = 0;

    while (file)
    {
        file.read(buffer.data(), buffer_size);
        std::streamsize count = file.gcount();
        if (count == 0)
        {
            break;
        }

        line_count
            += std::ranges::count(buffer.begin(), buffer.begin() + count, '\n');

        if (file.eof() && count > 0 && buffer[count - 1] != '\n')
        {
            line_count++;
        }
    }

    return line_count;
}

std::string parse_id(std::string_view line, char delimiter)
{
    auto parts = line | std::views::split(delimiter) | std::views::take(2);
    auto it = parts.begin();
    auto end = parts.end();

    if (it == end)
    {
        throw FileFormatException("failed to parse FID (empty line)");
    }
    std::string_view fid(*it);

    if (++it == end)
    {
        throw FileFormatException(
            "failed to parse FID and IID (missing delimiter)");
    }
    std::string_view iid(*it);
    return make_sample_id(fid, iid);
}

void parse_string(
    std::string_view line,
    std::vector<std::string_view>& out,
    size_t column_offset,
    char delimiter)
{
    out.clear();

    auto tokens
        = line | std::views::split(delimiter) | std::views::drop(column_offset);

    for (auto&& rng : tokens)
    {
        if (rng.empty())
        {
            throw DataParseException("empty value encountered");
        }
        out.emplace_back(rng);
    }
}

char detect_file_delimiter(std::ifstream& file)
{
    std::string probe_line;

    std::getline(file, probe_line);
    bool is_tab = !probe_line.empty() && probe_line.contains('\t');

    file.clear();
    file.seekg(0);

    return is_tab ? '\t' : ' ';
}

}  // namespace gelex::detail
