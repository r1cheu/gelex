#pragma once

#include <cstddef>
#include <span>
#include <string>

#include <fmt/color.h>
#include <fmt/ranges.h>

namespace gelex
{
std::string ToLowercase(std::string_view input);

std::string with_std(double value, double std);
std::string title(const std::string& text, size_t total_length = 80);
std::string subtitle(const std::string& str);

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
std::string
step(int n, int total, fmt::format_string<Args...> fmt_str, Args&&... args)
{
    return fmt::format("[{}/{}] ", n, total)
           + fmt::format(
               fmt::emphasis::bold, fmt_str, std::forward<Args>(args)...);
}

template <typename... Args>
std::string task(fmt::format_string<Args...> fmt_str, Args&&... args)
{
    return "  ● " + fmt::format(fmt_str, std::forward<Args>(args)...);
}

template <typename... Args>
std::string subtask(fmt::format_string<Args...> fmt_str, Args&&... args)
{
    return "    - " + fmt::format(fmt_str, std::forward<Args>(args)...);
}

template <typename... Args>
std::string success(fmt::format_string<Args...> fmt_str, Args&&... args)
{
    auto check_mark = fmt::format(fmt::fg(fmt::color::light_green), "  ✔ ");
    return check_mark + fmt::format(fmt_str, std::forward<Args>(args)...);
}

std::string format_names(
    std::span<const std::string> names,
    std::ptrdiff_t limit = 3);

/**
 * @brief Format the scale inverse chi-squared distribution
 * parameters into a human-readable string (e.g., "νS²χ⁻²(ν = nu, S² = s2)").
 *
 * @param nu Degrees of freedom.
 * @param s2 Scale parameter.
 * @return A formatted string representing the distribution parameters.
 */
std::string scale_inv_chisq(double nu, double s2);

std::string sigma_squared(const std::string& subscript);
std::string h2(const std::string& subscript);
std::string sigma_prior(const std::string& subscript, double nu, double s2);

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
