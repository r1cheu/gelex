#pragma once

#include <cstddef>
#include <span>
#include <string>
#include <vector>

#include <fmt/color.h>
#include <fmt/ranges.h>

namespace gelex
{
std::string ToLowercase(std::string_view input);

std::string with_std(double value, double std);
std::string title(const std::string& text, size_t total_length = 80);
std::string subtitle(const std::string& str);

std::string header_box(
    const std::string& title,
    const std::vector<std::pair<std::string, std::string>>& items,
    size_t width = 70);
std::string step_header(int current, int total, const std::string& description);
std::string separator(size_t width = 70, const std::string& c = "─");
std::string table_separator(size_t width = 70);

template <typename... Args>
std::string item(fmt::format_string<Args...> fmt_str, Args&&... args)
{
    return " ▪ " + fmt::format(fmt_str, std::forward<Args>(args)...);
}

template <typename... Args>
std::string subitem(fmt::format_string<Args...> fmt_str, Args&&... args)
{
    return "  - " + fmt::format(fmt_str, std::forward<Args>(args)...);
}

template <typename... Args>
std::string section(fmt::format_string<Args...> fmt_str, Args&&... args)
{
    return " "
           + fmt::format(
               fmt::emphasis::bold | fmt::fg(fmt::color::light_cyan),
               fmt_str,
               std::forward<Args>(args)...);
}

template <typename... Args>
std::string task(fmt::format_string<Args...> fmt_str, Args&&... args)
{
    return "   ● " + fmt::format(fmt_str, std::forward<Args>(args)...);
}

template <typename... Args>
std::string subtask(fmt::format_string<Args...> fmt_str, Args&&... args)
{
    return "    - " + fmt::format(fmt_str, std::forward<Args>(args)...);
}

template <typename... Args>
std::string subsubtask(fmt::format_string<Args...> fmt_str, Args&&... args)
{
    return "      └─ " + fmt::format(fmt_str, std::forward<Args>(args)...);
}

template <typename... Args>
std::string success(fmt::format_string<Args...> fmt_str, Args&&... args)
{
    auto check_mark = fmt::format(fmt::fg(fmt::color::light_green), "  ✔ ");
    return check_mark + fmt::format(fmt_str, std::forward<Args>(args)...);
}

template <typename... Args>
std::string warning_inline(fmt::format_string<Args...> fmt_str, Args&&... args)
{
    auto warning_mark
        = fmt::format(fmt::fg(fmt::color::orange), "  ! Warning: ");
    return warning_mark + fmt::format(fmt_str, std::forward<Args>(args)...);
}

template <typename... Args>
std::string progress_mark(fmt::format_string<Args...> fmt_str, Args&&... args)
{
    auto progress = fmt::format(fmt::fg(fmt::color::light_cyan), "  > ");
    return progress + fmt::format(fmt_str, std::forward<Args>(args)...);
}

std::string format_names(
    std::span<const std::string> names,
    std::ptrdiff_t limit = 3);

// fmt color might be removed
template <typename T>
auto green(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::rgb(0x00FF00)));
}

template <typename T>
auto pink(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::color::pink));
}

template <typename T>
auto cyan(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::color::cyan));
}

template <typename T>
auto rebecca_purple(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::color::rebecca_purple));
}
template <typename T>
auto rebecca_purple_vec(const T& value)
{
    return fmt::styled(
        fmt::join(value, ", "), fmt::fg(fmt::color::rebecca_purple));
}

template <typename T>
auto gray(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::color::gray));
}

template <typename T>
auto gold(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::color::gold));
}

template <typename T>
auto red(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::color::red));
}

template <typename T>
auto wine_red(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::rgb(0x982756)));
}
}  // namespace gelex
