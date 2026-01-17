#ifndef GELEX_UTILS_FORMATTER_H_
#define GELEX_UTILS_FORMATTER_H_

#include <cstddef>
#include <span>
#include <string>
#include <string_view>

#include <fmt/color.h>
#include <fmt/ranges.h>

namespace gelex
{

std::string header_box(
    std::string_view title,
    const std::span<std::pair<std::string, std::string>>& items,
    size_t width = 70);
std::string step_header(int current, int total, const std::string& description);
std::string separator(size_t width = 70, const std::string& c = "─");
std::string table_separator(size_t width = 70);

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
std::string progress_mark(fmt::format_string<Args...> fmt_str, Args&&... args)
{
    auto progress = fmt::format(fmt::fg(fmt::color::light_cyan), "  > ");
    return progress + fmt::format(fmt_str, std::forward<Args>(args)...);
}

std::string format_names(
    std::span<const std::string> names,
    std::ptrdiff_t limit = 3);

template <typename T>
auto rebecca_purple(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::color::rebecca_purple));
}
}  // namespace gelex

#endif  // GELEX_UTILS_FORMATTER_H_
