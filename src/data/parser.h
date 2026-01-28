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

#ifndef GELEX_DATA_PARSER_H_
#define GELEX_DATA_PARSER_H_

#include <concepts>
#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "gelex/exception.h"

namespace gelex::detail
{

template <typename T>
concept FileStream
    = std::derived_from<T, std::ios>
      && requires(std::filesystem::path p, std::ios::openmode m) {
             { T() };
         };

template <FileStream StreamType>
[[nodiscard]] StreamType open_file(
    const std::filesystem::path& path,
    std::ios::openmode mode,
    std::span<char> custom_buffer = {})
{
    if (std::filesystem::is_directory(path))
    {
        throw gelex::FileOpenException(
            std::format(
                "{}: is a directory, not a regular file", path.string()));
    }

    StreamType stream;
    if (!custom_buffer.empty())
    {
        stream.rdbuf()->pubsetbuf(
            custom_buffer.data(),
            static_cast<std::streamsize>(custom_buffer.size()));
    }

    stream.open(path, mode);

    if (!stream.is_open())
    {
        if ((mode & std::ios::in) && !std::filesystem::exists(path))
        {
            throw FileNotFoundException(
                std::format("{}: not found", path.string()));
        }
        throw FileOpenException(
            std::format("{}: failed to open file", path.string()));
    }

    if ((mode & std::ios::in) && std::filesystem::is_regular_file(path))
    {
        std::error_code ec;
        if (std::filesystem::file_size(path, ec) == 0 && !ec)
        {
            throw FileFormatException(
                std::format("{}: is empty", path.string()));
        }
    }

    return stream;
}

size_t estimate_line_count(
    const std::filesystem::path& path,
    int sample_lines = 5);

// counts total number of lines in a file
size_t count_total_lines(const std::filesystem::path& path);

size_t count_num_columns(std::string_view line, char delimiter = '\t');

template <typename T = double>
T parse_number(std::string_view sv)
{
    if (sv.empty())
    {
        throw NumberParseException(
            std::format("empty string cannot be parsed as number"));
    }

    T value{};
    const char* first = sv.data();
    const char* last = sv.data() + sv.size();
    auto [ptr, ec] = std::from_chars(first, last, value);

    if (ec == std::errc() && ptr == last)
    {
        return value;
    }
    throw NumberParseException(
        std::format("failed to parse '{}' as number", sv));
}

double parse_nth_double(
    std::string_view line,
    size_t column_index,
    char delimiter = '\t');

void parse_all_doubles(
    std::string_view line,
    std::vector<double>& out,
    size_t column_offset = 0,
    char delimiter = '\t');

std::vector<std::string_view> parse_header(
    std::string_view line,
    char delimiter = '\t');

std::string
parse_id(std::string_view line, bool iid_only, char delimiter = '\t');

char detect_file_delimiter(std::ifstream& file);

void parse_string(
    std::string_view line,
    std::vector<std::string_view>& out,
    size_t column_offset = 0,
    char delimiter = '\t');

}  // namespace gelex::detail

#endif  // GELEX_DATA_PARSER_H_
