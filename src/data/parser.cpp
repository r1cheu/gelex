#include "parser.h"

#include <expected>
#include <ranges>
#include <string_view>
#include <vector>

#include "gelex/error.h"

namespace gelex
{
namespace detail
{

auto parse_nth_double(
    std::string_view line,
    size_t column_index,
    std::string_view delimiters) noexcept -> std::expected<double, Error>
{
    auto token_view = line | std::views::split(delimiters)
                      | std::views::filter([](auto&& r) { return !r.empty(); })
                      | std::views::drop(column_index) | std::views::take(1);

    if (auto it = token_view.begin(); it != token_view.end())
    {
        return try_parse_double(*it);
    }
    return std::unexpected(
        Error{
            .code = ErrorCode::InvalidRange,
            .message
            = std::format("column index {} is out of range", column_index)});
}

auto parse_id(
    std::string_view line,
    bool iid_only,
    std::string_view delimiters) noexcept -> std::expected<std::string, Error>
{
    auto ids = line | std::views::split(delimiters) | std::views::take(2)
               | std::views::filter([](auto&& sr) { return !sr.empty(); })
               | std::views::transform([](auto&& subrange)
                                       { return std::string_view(subrange); })
               | std::ranges::to<std::vector<std::string_view>>();

    if (ids.size() < 2)
    {
        return std::unexpected(
            Error{
                .code = ErrorCode::InvalidFile,
                .message = "failed to parse FID and IID"});
    }

    if (iid_only)
    {
        return std::string(ids[1]);
    }

    return std::format("{}_{}", ids[0], ids[1]);
}

size_t count_num_columns(
    std::string_view line,
    std::string_view delimiters) noexcept
{
    auto tokens = line | std::views::split(delimiters)
                  | std::views::filter([](auto&& subrange)
                                       { return !subrange.empty(); });

    return static_cast<size_t>(std::ranges::distance(tokens));
}

auto parse_header(std::string_view line, std::string_view delimiters) noexcept
    -> std::expected<std::vector<std::string_view>, Error>
{
    auto header = detail::parse_string(line, 0, delimiters);
    if (header.empty())
    {
        return std::unexpected(
            Error{.code = ErrorCode::WrongHeader, .message = "empty header"});
    }

    if (header.size() < 2)
    {
        return std::unexpected(
            Error{
                .code = ErrorCode::WrongHeader,
                .message = "header contains less than 2 columns."});
    }

    if (header[0] != "FID" || header[1] != "IID")
    {
        return std::unexpected(
            Error{
                .code = ErrorCode::WrongHeader,
                .message = std::format(
                    "first two columns are {} and {}, but "
                    "expected 'FID' and 'IID'.",
                    header[0],
                    header[1])});
    }
    return header;
}

auto parse_string(
    std::string_view line,
    size_t column_offset,
    std::string_view delimiters) noexcept -> std::vector<std::string_view>
{
    auto tokens
        = line | std::views::split(delimiters)
          | std::views::filter([](auto&& sr) { return !sr.empty(); })
          | std::views::drop(column_offset)
          | std::views::transform([](auto&& subrange)
                                  { return std::string_view(subrange); });

    return std::ranges::to<std::vector<std::string_view>>(tokens);
}

auto parse_all_doubles(
    std::string_view line,
    size_t column_offset,
    std::string_view delimiters) noexcept
    -> std::expected<std::vector<double>, Error>
{
    auto token_view = line | std::views::split(delimiters)
                      | std::views::filter([](auto&& r) { return !r.empty(); })
                      | std::views::drop(column_offset);

    std::vector<double> result;
    int column_index{static_cast<int>(column_offset)};
    for (auto&& token_range : token_view)
    {
        auto parsed_value = try_parse_double(token_range);
        if (parsed_value)
        {
            result.push_back(*parsed_value);
        }
        else
        {
            auto& error = parsed_value.error();
            error.message
                = std::format("{} at column {}", error.message, column_index);
            return std::unexpected(parsed_value.error());
        }
        column_index++;
    }
    return result;
}

}  // namespace detail
}  // namespace gelex
