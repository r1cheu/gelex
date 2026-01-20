#ifndef GELEX_UTILS_FORMATTER_H_
#define GELEX_UTILS_FORMATTER_H_

#include <Eigen/Core>
#include <cmath>
#include <cstddef>
#include <span>
#include <string>
#include <string_view>

#include <fmt/base.h>
#include <fmt/color.h>
#include <fmt/ranges.h>

#include "../src/types/freq_effect.h"

namespace gelex
{

std::string header_box(
    std::string_view title,
    const std::span<std::pair<std::string, std::string>>& items,
    size_t width = 70);
std::string step_header(int current, int total, const std::string& description);
std::string separator(size_t width = 70, const std::string& c = "─");
std::string table_separator(size_t width = 70);

std::string format_eta(double seconds);

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

struct HumanReadable
{
    explicit HumanReadable(size_t value) : value(static_cast<double>(value)) {}
    explicit HumanReadable(Eigen::Index value)
        : value(static_cast<double>(value))
    {
    }
    explicit HumanReadable(double value) : value(value) {}
    double value;
};

}  // namespace gelex

namespace fmt
{
template <>
struct formatter<gelex::freq::GrmType> : formatter<string_view>
{
    // parse is inherited from formatter<string_view>.

    auto format(gelex::freq::GrmType t, format_context& ctx) const
        -> format_context::iterator;
};

template <>
struct formatter<gelex::HumanReadable> : formatter<double>
{
    auto format(gelex::HumanReadable hr, format_context& ctx) const
        -> format_context::iterator;
};
}  // namespace fmt

#endif  // GELEX_UTILS_FORMATTER_H_
