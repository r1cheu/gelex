#include "formatter.h"

#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include <fmt/chrono.h>
#include <fmt/format.h>

namespace gelex
{
std::string with_std(double value, double std)
{
    return fmt::format("{:.6f} \u00B1 {:.4f}", value, std);
}

std::string title(const std::string& text, size_t total_length)
{
    if (text.empty())
    {
        return std::string(total_length, '=');  // Fallback to ASCII hyphen
    }

    size_t text_length = text.length();
    if (text_length >= total_length)
    {
        return text;
    }

    size_t remaining_space = total_length - text_length;
    size_t left_hyphens = remaining_space / 2;
    size_t right_hyphens = remaining_space - left_hyphens;  // Handles odd cases

    std::string left(left_hyphens, '=');
    std::string right(right_hyphens, '=');
    return left + text + right;
}

std::string format_names(
    std::span<const std::string> names,
    std::ptrdiff_t limit)
{
    if (names.empty())
    {
        return "None";
    }
    if (names.size() <= static_cast<size_t>(limit))
    {
        return fmt::format("({})", fmt::join(names, ", "));
    }

    std::span<const std::string> sub = names.first(limit - 1);
    return fmt::format(
        "({}, ... and {} more)",
        fmt::join(sub, ", "),
        names.size() - (limit - 1));
}

std::string subtitle(const std::string& str)
{
    return fmt::format("{}", wine_red("[" + str + "]"));
}

std::string ToLowercase(std::string_view input)
{
    std::string result(input.begin(), input.end());
    std::transform(result.begin(), result.end(), result.begin(), ::tolower);
    return result;
}
}  // namespace gelex
