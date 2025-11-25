#pragma once

#include <charconv>
#include <concepts>
#include <filesystem>
#include <format>
#include <fstream>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

#include "gelex/exception.h"

namespace gelex::detail
{

template <typename T>
concept FileStream
    = std::derived_from<T, std::ios_base>
      && requires(std::filesystem::path p, std::ios_base::openmode m) {
             { T() };
         };

template <FileStream StreamType>
[[nodiscard]] StreamType open_file(
    const std::filesystem::path& path,
    std::ios_base::openmode mode,
    std::span<char> custom_buffer = {})
{
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
        if ((mode & std::ios_base::in) && !std::filesystem::exists(path))
        {
            throw FileNotFoundException(path);
        }

        throw FileIOException(
            enrich_with_file_info("Failed to open file", path));
    }

    if ((mode & std::ios_base::in) && std::filesystem::is_regular_file(path)
        && std::filesystem::file_size(path) == 0)
    {
        throw InvalidFileException(
            enrich_with_file_info("File is empty", path));
    }

    return stream;
}

template <std::ranges::contiguous_range R>
double try_parse_double(R&& token_range)
{
    const char* begin = std::to_address(std::ranges::begin(token_range));
    const char* end = std::to_address(std::ranges::end(token_range));

    double value{};
    const auto result = std::from_chars(begin, end, value);

    if (result.ec == std::errc() && result.ptr == end)
    {
        return value;
    }

    throw NotNumberException(
        std::format(
            "failed to parse '{}' as a number",
            std::string_view(begin, end - begin)));
}

size_t count_total_line(const std::filesystem::path& path);

size_t count_total_lines(const std::filesystem::path& path);

double parse_nth_double(
    std::string_view line,
    size_t column_index,
    char delimiter = '\t');

std::string
parse_id(std::string_view line, bool iid_only, char delimiter = '\t');

std::vector<std::string_view> parse_header(
    std::string_view line,
    char delimiter = '\t');

size_t count_num_columns(std::string_view line, char delimiter = '\t');

std::vector<std::string_view> parse_string(
    std::string_view line,
    size_t column_offset = 0,
    char delimiter = '\t');

std::vector<double> parse_all_doubles(
    std::string_view line,
    size_t column_offset = 0,
    char delimiter = '\t');

}  // namespace gelex::detail
