#include "parser.h"

#include <algorithm>
#include <cstddef>

namespace gelex::detail
{

size_t count_total_lines(const std::filesystem::path& path)
{
    auto file = open_file<std::ifstream>(
        path, std::ios_base::in | std::ios_base::binary);

    constexpr auto buffer_size = static_cast<const size_t>(64 * 1024);
    std::vector<char> buffer(buffer_size);

    size_t line_count = 0;
    while (file.read(buffer.data(), buffer.size()))
    {
        line_count += std::ranges::count(buffer, '\n');
    }

    if (file.gcount() > 0)
    {
        std::span<char> last_chunk(buffer.data(), file.gcount());
        line_count += std::ranges::count(last_chunk, '\n');

        if (last_chunk.back() != '\n')
        {
            line_count++;
        }
    }

    return line_count;
}

auto split_line(std::string_view line, char delimiter)
{
    return line | std::views::split(delimiter)
           | std::views::filter([](auto&& r)
                                { return !std::ranges::empty(r); });
}

double
parse_nth_double(std::string_view line, size_t column_index, char delimiter)
{
    auto tokens = split_line(line, delimiter) | std::views::drop(column_index)
                  | std::views::take(1);

    if (auto it = tokens.begin(); it != tokens.end())
    {
        return try_parse_double(*it);
    }

    throw InvalidRangeException(
        std::format("column index {} is out of range", column_index));
}

std::string parse_id(std::string_view line, bool iid_only, char delimiter)
{
    auto tokens = split_line(line, delimiter) | std::views::take(2);

    std::vector<std::string_view> ids;
    ids.reserve(2);

    for (auto&& sub : tokens)
    {
        ids.emplace_back(sub.begin(), sub.end());
    }

    if (ids.size() < 2)
    {
        throw InvalidFileException("failed to parse FID and IID");
    }

    if (iid_only)
    {
        return std::string(ids[1]);
    }

    return std::format("{}_{}", ids[0], ids[1]);
}

size_t count_num_columns(std::string_view line, char delimiter)
{
    auto tokens = split_line(line, delimiter);
    return static_cast<size_t>(std::ranges::distance(tokens));
}

std::vector<std::string_view>
parse_string(std::string_view line, size_t column_offset, char delimiter)
{
    auto tokens_view
        = split_line(line, delimiter) | std::views::drop(column_offset);

    std::vector<std::string_view> result;
    for (auto&& rng : tokens_view)
    {
        result.emplace_back(rng.begin(), rng.end());
    }
    return result;
}

std::vector<std::string_view> parse_header(
    std::string_view line,
    char delimiter)
{
    auto header = parse_string(line, 0, delimiter);

    if (header.empty())
    {
        throw WrongHeaderException("empty header");
    }

    if (header.size() < 2)
    {
        throw WrongHeaderException("header contains less than 2 columns.");
    }

    if (header[0] != "FID" || header[1] != "IID")
    {
        throw WrongHeaderException(
            std::format(
                "first two columns are '{}' and '{}', but expected 'FID' and "
                "'IID'.",
                header[0],
                header[1]));
    }
    return header;
}

std::vector<double>
parse_all_doubles(std::string_view line, size_t column_offset, char delimiter)
{
    auto tokens_view
        = split_line(line, delimiter) | std::views::drop(column_offset);

    std::vector<double> result;
    int column_index = static_cast<int>(column_offset);
    for (auto&& token_range : tokens_view)
    {
        try
        {
            result.push_back(try_parse_double(token_range));
        }
        catch (const NotNumberException& ex)
        {
            throw InvalidDataException(
                std::format("{} at column {}", ex.what(), column_index));
        }
        column_index++;
    }
    return result;
}

}  // namespace gelex::detail
