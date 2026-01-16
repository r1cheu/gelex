#include "formatter.h"

#include <algorithm>
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

std::string header_box(
    const std::string& title,
    const std::vector<std::pair<std::string, std::string>>& items,
    size_t width)
{
    std::string result;
    const std::string h_line = "─";
    const std::string v_line = "│";

    // Top border
    result += fmt::format(fmt::fg(fmt::color::light_cyan), "┌── {} ", title);
    size_t title_len = 4 + title.length() + 1;  // "┌── " + title + " "
    if (title_len < width - 1)
    {
        for (size_t i = 0; i < width - title_len - 1; ++i)
        {
            result
                += fmt::format(fmt::fg(fmt::color::light_cyan), "{}", h_line);
        }
        result += fmt::format(fmt::fg(fmt::color::light_cyan), "┐\n");
    }
    else
    {
        result += fmt::format(fmt::fg(fmt::color::light_cyan), "┐\n");
    }

    // Items
    for (const auto& [key, value] : items)
    {
        std::string line = fmt::format(
            "  {}{}  : {}", key, std::string(12 - key.length(), ' '), value);
        size_t line_len = 2 + key.length() + (12 - key.length()) + 5
                          + value.length();  // visible length
        if (line_len < width - 2)
        {
            result += fmt::format(
                fmt::fg(fmt::color::light_cyan),
                "{}{}{}│\n",
                v_line,
                line,
                std::string(width - line_len - 2, ' '));
        }
        else
        {
            result += fmt::format(
                fmt::fg(fmt::color::light_cyan), "{}{}│\n", v_line, line);
        }
    }

    // Bottom border
    result += fmt::format(fmt::fg(fmt::color::light_cyan), "└");
    for (size_t i = 0; i < width - 2; ++i)
    {
        result += fmt::format(fmt::fg(fmt::color::light_cyan), "{}", h_line);
    }
    result += fmt::format(fmt::fg(fmt::color::light_cyan), "┘");

    return result;
}

std::string step_header(int current, int total, const std::string& description)
{
    return fmt::format(
        fmt::emphasis::bold | fmt::fg(fmt::color::light_cyan),
        "[{}/{}] {}",
        current,
        total,
        description);
}

std::string separator(size_t width, const std::string& c)
{
    std::string result;
    for (size_t i = 0; i < width; ++i)
    {
        result += c;
    }
    return result;
}

std::string table_separator(size_t width)
{
    return "  " + separator(width - 2);
}
}  // namespace gelex
