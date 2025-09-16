#pragma once

#include <charconv>
#include <expected>
#include <filesystem>
#include <format>
#include <ranges>
#include <string_view>
#include <vector>

namespace gelex
{
namespace detail
{

enum class ParseError : uint8_t
{
    NotNumber,
    InvalidColumn,
    InvalidFile,
};

template <std::ranges::range R>
auto try_parse_double(const R& token_range) noexcept
    -> std::expected<double, ParseError>
{
    const std::string_view token{token_range};
    double value{};
    const auto result
        = std::from_chars(token.data(), token.data() + token.size(), value);

    if (result.ec == std::errc() && result.ptr == token.data() + token.size())
    {
        return value;
    }
    return std::unexpected(ParseError::NotNumber);
}

}  // namespace detail
}  // namespace gelex

// namespace std
namespace std
{
template <>
struct formatter<gelex::detail::ParseError> : std::formatter<std::string_view>
{
    auto format(gelex::detail::ParseError p, format_context& ctx) const
    {
        std::string_view name = "Unknown Error";
        switch (p)
        {
            case gelex::detail::ParseError::NotNumber:
                name = "NotNumber";
                break;
            case gelex::detail::ParseError::InvalidColumn:
                name = "InvalidColumn";
                break;
            case gelex::detail::ParseError::InvalidFile:
                name = "InvalidFile";
                break;
        }
        return formatter<string_view>::format(name, ctx);
    }
};
}  // namespace std

namespace gelex
{
namespace detail
{

/**
 * @brief Parse the nth double from a delimited string.
 *
 * @param line The input string to parse
 * @param column_index column index (0-based)
 * @param delimiters
 */
auto parse_nth_double(
    std::string_view line,
    size_t column_index,
    std::string_view delimiters = "\t") noexcept
    -> std::expected<double, ParseError>;

auto parse_id(
    std::string_view line,
    bool iid_only,
    std::string_view delimiters = "\t") noexcept
    -> std::expected<std::string, ParseError>;

auto parse_header(
    std::string_view line,
    const std::filesystem::path& path,
    std::string_view delimiters = "\t") -> std::vector<std::string_view>;

/**
 * @brief Parses a string into tokens based on specified delimiters. Return
 * vector<string_view> for speed, should be used carefully.
 *
 * @param line The input string to parse.
 * @param column_offset The starting position in the string to begin
 * parsing.
 * @param delimiters A string containing delimiter characters to split the
 * input.
 * @return std::vector<std::string_view> A vector of string views
 * representing the parsed tokens.
 */
auto parse_string(
    std::string_view line,
    size_t column_offset = 0,
    std::string_view delimiters = "\t") noexcept -> std::vector<std::string>;

/**
 * @brief Parses all double values from a string, starting at a given column
 * offset and using specified delimiters.
 *
 * This function splits the input string into tokens using the provided
 * delimiters, starting from the specified column offset, and attempts to parse
 * each token as a double. If all tokens are successfully parsed, returns a
 * vector of doubles. If any token cannot be parsed as a double, returns a
 * ParseError.
 *
 * @param line The input string to parse.
 * @param column_offset The starting position in the string to begin parsing.
 * @param delimiters A string containing delimiter characters to split the
 * input.
 */
auto parse_all_doubles(
    std::string_view line,
    size_t column_offset = 0,
    std::string_view delimiters = "\t") noexcept
    -> std::expected<std::vector<double>, ParseError>;
}  // namespace detail
}  // namespace gelex
