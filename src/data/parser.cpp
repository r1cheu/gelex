#include "parser.h"

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <ranges>
#include <vector>

namespace gelex::detail
{

size_t count_total_lines(const std::filesystem::path& path)
{
    constexpr size_t buffer_size = 128 * 1024;

    std::vector<char> buffer(buffer_size);

    auto file = open_file<std::ifstream>(
        path, std::ios_base::in | std::ios_base::binary);

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

std::string_view get_nth_token(std::string_view line, size_t n, char delimiter)
{
    auto view = line | std::views::split(delimiter) | std::views::drop(n);
    if (view.begin() == view.end())
    {
        return {};
    }
    return std::string_view(*view.begin());
}

double
parse_nth_double(std::string_view line, size_t column_index, char delimiter)
{
    std::string_view token = get_nth_token(line, column_index, delimiter);

    if (!token.empty() || (line.size() > 0 && column_index == 0))
    {
        return parse_number<double>(token);
    }

    throw ColumnRangeException(
        std::format("Column {} is out of range", column_index));
}

std::string parse_id(std::string_view line, bool iid_only, char delimiter)
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
    if (iid_only)
    {
        return std::string(iid);
    }
    return std::format("{}_{}", fid, iid);
}

size_t count_num_columns(std::string_view line, char delimiter)
{
    if (line.empty())
    {
        return 0;
    }
    return static_cast<size_t>(std::ranges::count(line, delimiter)) + 1;
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

std::vector<std::string_view> parse_header(
    std::string_view line,
    char delimiter)
{
    std::vector<std::string_view> header;
    header.reserve(16);
    parse_string(line, header, 0, delimiter);

    if (header.size() < 2)
    {
        throw HeaderFormatException(
            std::format("header contains only {} columns.", header.size()));
    }

    if (header[0] != "FID" || header[1] != "IID")
    {
        throw HeaderFormatException(
            std::format(
                "first two columns are '{}' and '{}', expected 'FID' and "
                "'IID'.",
                header[0],
                header[1]));
    }
    return header;
}

void parse_all_doubles(
    std::string_view line,
    std::vector<double>& out,
    size_t column_offset,
    char delimiter)
{
    out.clear();
    auto tokens = line | std::views::split(delimiter)
                  | std::views::drop(column_offset) | std::views::enumerate;

    for (auto&& [idx, rng] : tokens)
    {
        std::string_view token(rng);

        if (token.empty())
        {
            continue;
        }

        try
        {
            out.push_back(parse_number<double>(token));
        }
        catch (const NumberParseException& ex)
        {
            size_t actual_col = idx + column_offset;
            throw DataParseException(
                std::format("{} at column {}", ex.what(), actual_col));
        }
    }
}

}  // namespace gelex::detail
