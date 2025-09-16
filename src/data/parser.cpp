#include "parser.h"

#include <expected>
#include <ranges>
#include <string_view>
#include <vector>

namespace gelex
{
namespace detail
{

auto parse_nth_double(
    std::string_view line,
    size_t column_index,
    std::string_view delimiters) noexcept -> std::expected<double, ParseError>
{
    auto token_view = line | std::views::split(delimiters)
                      | std::views::filter([](auto&& r) { return !r.empty(); })
                      | std::views::drop(column_index);

    if (auto it = token_view.begin(); it != token_view.end())
    {
        return try_parse_double(*it);
    }
    return std::unexpected(ParseError::InvalidColumn);
}

auto parse_id(
    std::string_view line,
    bool iid_only,
    std::string_view delimiters) noexcept
    -> std::expected<std::string, ParseError>
{
    auto ids = line | std::views::split(delimiters) | std::views::take(2)
               | std::views::filter([](auto&& sr) { return !sr.empty(); })
               | std::views::transform([](auto&& subrange)
                                       { return std::string_view(subrange); })
               | std::ranges::to<std::vector<std::string_view>>();

    if (ids.size() < 2)
    {
        return std::unexpected(ParseError::InvalidFile);
    }

    if (iid_only)
    {
        return std::string(ids[1]);
    }

    return std::format("{}_{}", ids[0], ids[1]);
}

std::vector<std::string_view> parse_header(
    std::string_view line,
    const std::filesystem::path& path,
    std::string_view delimiters)
{
    auto tokens = detail::parse_string(line, 0, delimiters);
    if (tokens.empty())
    {
        if (!line.empty())
        {
            throw std::runtime_error(
                std::format(
                    "Invalid header in file [{}]: Header is empty or contains "
                    "only whitespace.",
                    path.string()));
        }
        throw std::runtime_error(
            std::format(
                "Invalid header in file [{}]: Header is empty.",
                path.string()));
    }

    if (tokens.size() < 2)
    {
        throw std::runtime_error(
            std::format(
                "Invalid header in file [{}]: Header contains fewer than two "
                "columns.",
                path.string()));
    }

    if (tokens[0] != "FID" || tokens[1] != "IID")
    {
        throw std::runtime_error(
            std::format(
                "Invalid header in file [{}]: First two columns must be 'FID' "
                "and 'IID', but found '{}' and '{}'.",
                path.string(),
                tokens[0],
                tokens[1]));
    }

    return tokens;
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
    -> std::expected<std::vector<double>, ParseError>
{
    auto token_view = line | std::views::split(delimiters)
                      | std::views::filter([](auto&& r) { return !r.empty(); })
                      | std::views::drop(column_offset);

    std::vector<double> result;
    for (auto&& token_range : token_view)
    {
        auto parsed_value = try_parse_double(token_range);
        if (parsed_value)
        {
            result.push_back(*parsed_value);
        }
        else
        {
            return std::unexpected(parsed_value.error());
        }
    }
    return result;
}

}  // namespace detail
}  // namespace gelex
