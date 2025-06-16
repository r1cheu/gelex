#include "gelex/utils/formatter.h"

#include <string>
#include <string_view>

#include <fmt/chrono.h>
#include <fmt/format.h>
#include <armadillo>

namespace gelex
{
std::string format_sigma_squared(const std::string& name)
{
    return fmt::format("\u03C3\u00B2_{}", name);
}

std::string format_value_with_std(double value, double std)
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

std::string subtitle(const std::string& str)
{
    return fmt::format("{}", wine_red("[" + str + "]"));
}

std::string item(const std::string& item)
{
    return fmt::format(" \u25AA {}", item);
}

std::string subitem(const std::string& item)
{
    return fmt::format(" - {}", item);
}

std::string ToLowercase(std::string_view input)
{
    std::string result(input.begin(), input.end());
    std::transform(result.begin(), result.end(), result.begin(), ::tolower);
    return result;
}
}  // namespace gelex
