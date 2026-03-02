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

#include "gelex/infra/utils/formatter.h"

#include <fmt/chrono.h>
#include <fmt/format.h>
#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace gelex
{

std::vector<std::string> header_box(
    const std::string_view title,
    const std::span<std::pair<std::string, std::string>>& items,
    size_t width)
{
    std::vector<std::string> lines;
    const std::string h_line = "─";
    const std::string v_line = "│";

    lines.emplace_back("");

    // Top border
    std::string top
        = fmt::format(fmt::fg(fmt::color::light_cyan), "┌── {} ", title);
    size_t title_len = 4 + title.length() + 1;  // "┌── " + title + " "
    if (title_len < width - 1)
    {
        for (size_t i = 0; i < width - title_len - 1; ++i)
        {
            top += fmt::format(fmt::fg(fmt::color::light_cyan), "{}", h_line);
        }
    }
    top += fmt::format(fmt::fg(fmt::color::light_cyan), "┐");
    lines.push_back(std::move(top));

    // Items
    for (const auto& [key, value] : items)
    {
        std::string content = fmt::format(
            "  {}{}  : {}", key, std::string(12 - key.length(), ' '), value);
        size_t line_len = 2 + key.length() + (12 - key.length()) + 5
                          + value.length();  // visible length
        std::string row;
        if (line_len < width - 2)
        {
            row = fmt::format(
                fmt::fg(fmt::color::light_cyan),
                "{}{}{}│",
                v_line,
                content,
                std::string(width - line_len - 1, ' '));
        }
        else
        {
            row = fmt::format(
                fmt::fg(fmt::color::light_cyan), "{}{}│", v_line, content);
        }
        lines.push_back(std::move(row));
    }

    // Bottom border
    std::string bottom = fmt::format(fmt::fg(fmt::color::light_cyan), "└");
    for (size_t i = 0; i < width - 2; ++i)
    {
        bottom += fmt::format(fmt::fg(fmt::color::light_cyan), "{}", h_line);
    }
    bottom += fmt::format(fmt::fg(fmt::color::light_cyan), "┘");
    lines.push_back(std::move(bottom));

    return lines;
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

std::string named_section(std::string_view name, size_t width, size_t indent)
{
    std::string result(indent, ' ');
    result += "── ";
    result += name;
    result += " ";
    size_t used = indent + 3 + name.size() + 1;
    size_t remaining = used < width ? width - used : 0;
    result += separator(remaining);
    return fmt::format(fmt::emphasis::bold, "{}", result);
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
auto fmt::formatter<gelex::AbbrNumber>::format(
    gelex::AbbrNumber hr,
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
