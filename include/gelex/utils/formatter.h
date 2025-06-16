#pragma once

#include <string>

#include <fmt/color.h>
#include <fmt/ranges.h>

namespace gelex
{
std::string ToLowercase(std::string_view input);

// help functions for logging
std::string format_sigma_squared(const std::string& name);
std::string format_value_with_std(double value, double std);
std::string title(const std::string& text, size_t total_length = 80);
std::string subtitle(const std::string& str);
std::string item(const std::string& item);
std::string subitem(const std::string& item);

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
