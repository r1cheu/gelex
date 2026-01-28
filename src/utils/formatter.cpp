/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "formatter.h"

#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include "../src/types/freq_effect.h"

#include <fmt/chrono.h>
#include <fmt/format.h>

namespace gelex
{

std::string header_box(
    const std::string_view title,
    const std::span<std::pair<std::string, std::string>>& items,
    size_t width)
{
    std::string result = "\n";
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
                std::string(width - line_len - 1, ' '));
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

std::string format_eta(double seconds)
{
    if (seconds < 0 || seconds > 3600 * 60 * 99)
    {
        return "--:--:--";
    }
    int total_seconds = static_cast<int>(seconds);
    int h = total_seconds / 3600;
    int m = total_seconds / 60;
    int s = total_seconds % 60;
    return std::format("{:02d}:{:02d}:{:02d}", h, m, s);
}

std::string format_names(
    std::span<const std::string> names,
    std::ptrdiff_t limit)
{
    if (names.empty())
    {
        return "";
    }
    if (names.size() <= static_cast<size_t>(limit))
    {
        return fmt::format("{}", fmt::join(names, ", "));
    }
    std::span<const std::string> displayed = names.first(limit);
    std::ptrdiff_t remaining
        = static_cast<std::ptrdiff_t>(names.size()) - limit;
    return fmt::format(
        "{}, ... and {} more", fmt::join(displayed, ", "), remaining);
}
}  // namespace gelex

namespace fmt
{
auto formatter<gelex::freq::GrmType>::format(
    gelex::freq::GrmType t,
    format_context& ctx) const -> format_context::iterator
{
    string_view name = "unknown";
    switch (t)
    {
        case (gelex::freq::GrmType::A):
            name = "A";
            break;
        case (gelex::freq::GrmType::D):
            name = "D";
            break;
        case (gelex::freq::GrmType::AA):
            name = "AA";
            break;
        case (gelex::freq::GrmType::AD):
            name = "AD";
            break;
        case (gelex::freq::GrmType::DD):
            name = "DD";
            break;
        case (gelex::freq::GrmType::Unknown):
        default:
            break;
    }
    return formatter<string_view>::format(name, ctx);
}

auto fmt::formatter<gelex::HumanReadable>::format(
    gelex::HumanReadable hr,
    format_context& ctx) const -> format_context::iterator
{
    double v = hr.value;
    const char* suffix = "";

    if (v >= 1e9)
    {
        v /= 1e9;
        suffix = "G";
    }
    else if (v >= 1e6)
    {
        v /= 1e6;
        suffix = "M";
    }
    else if (v >= 1e3)
    {
        v /= 1e3;
        suffix = "k";
    }

    if (v >= 100 || std::abs(v - std::round(v)) < 0.05)
    {
        return fmt::format_to(ctx.out(), "{:.0f}{}", v, suffix);
    }
    return fmt::format_to(ctx.out(), "{:.1f}{}", v, suffix);
}

}  // namespace fmt
