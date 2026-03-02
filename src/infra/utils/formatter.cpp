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

namespace gelex
{

std::string
command_banner(std::string_view version, std::string_view task, size_t width)
{
    // ▐ GELEX ▌  (9 visible chars, bold + light_cyan)
    std::string badge = fmt::format(
        fmt::emphasis::bold | fmt::fg(fmt::color::light_cyan),
        "\u258c GELEX \u2590");

    // "  v{version}  " plain text; visible = version.size() + 5
    std::string ver_str = fmt::format("  v{}  ", version);

    // ──  {task} {fill}  (bold); prefix visible = 4 + task.size() + 1
    size_t visible = 9 + (version.size() + 5) + 4 + task.size() + 1;
    size_t fill_count = visible < width ? width - visible - 1 : 0;
    std::string fill;
    fill.reserve(fill_count * 3);  // U+2500 encodes as 3 UTF-8 bytes
    for (size_t i = 0; i < fill_count; ++i)
    {
        fill += "\u2500";
    }
    std::string sep_task = fmt::format(
        fmt::emphasis::bold, "\u2500\u2500  {} {}\u256e", task, fill);

    return badge + ver_str + sep_task;
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

std::string done_message(double elapsed_seconds)
{
    auto sparkle = fmt::format(fmt::fg(fmt::color::yellow), "✨");
    return fmt::format("╰── {} Done in {:.2f}s", sparkle, elapsed_seconds);
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
