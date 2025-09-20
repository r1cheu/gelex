#pragma once

#include <memory>
#include <string_view>

#include <fmt/color.h>
#include <fmt/ranges.h>
#include <spdlog/formatter.h>
#include <spdlog/logger.h>

namespace gelex
{
namespace detail
{

std::string ToLowercase(std::string_view input);

std::string with_std(double value, double std);
std::string title(const std::string& text, size_t total_length = 80);
std::string subtitle(const std::string& str);

template <typename... Args>
std::string item(fmt::format_string<Args...> fmt_str, Args&&... args)
{
    return " \u25AA " + fmt::format(fmt_str, std::forward<Args>(args)...);
}

template <typename... Args>
std::string subitem(fmt::format_string<Args...> fmt_str, Args&&... args)
{
    return "  - " + fmt::format(fmt_str, std::forward<Args>(args)...);
}
std::string sigma_prior(std::string_view subscript, double nu, double s2);
std::string sigma_squared(std::string_view subscript);
std::string h2(std::string_view subscript);

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
}  // namespace detail

}  // namespace gelex
