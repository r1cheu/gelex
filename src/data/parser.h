#pragma once

#include <charconv>
#include <expected>
#include <format>
#include <ranges>
#include <string_view>
#include <vector>

#include "gelex/error.h"

namespace gelex
{
namespace detail
{

template <std::ranges::range R>
auto try_parse_double(const R& token_range) noexcept
    -> std::expected<double, Error>
{
    const std::string_view token{token_range};
    double value{};
    const auto result
        = std::from_chars(token.data(), token.data() + token.size(), value);

    if (result.ec == std::errc() && result.ptr == token.data() + token.size())
    {
        return value;
    }
    return std::unexpected(
        Error{
            .code = ErrorCode::NotNumber,
            .message = std::format("failed to parse '{}' as a number", token)});
}

}  // namespace detail
}  // namespace gelex

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
    -> std::expected<double, Error>;

auto parse_id(
    std::string_view line,
    bool iid_only,
    std::string_view delimiters = "\t") noexcept
    -> std::expected<std::string, Error>;

auto parse_header(
    std::string_view line,
    std::string_view delimiters = "\t") noexcept
    -> std::expected<std::vector<std::string_view>, Error>;

size_t count_num_columns(
    std::string_view line,
    std::string_view delimiters = "\t") noexcept;

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
    std::string_view delimiters = "\t") noexcept
    -> std::vector<std::string_view>;

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
    -> std::expected<std::vector<double>, Error>;
}  // namespace detail
}  // namespace gelex
