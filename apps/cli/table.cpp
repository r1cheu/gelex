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

#include "cli/table.h"

#include <algorithm>
#include <cstddef>
#include <fmt/color.h>
#include <fmt/format.h>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "cli/formatter.h"

namespace cli
{

namespace
{

constexpr std::string_view indent = "   ";
constexpr std::string_view gap = "   ";

auto display_width(std::string_view s) -> size_t
{
    size_t width = 0;
    for (size_t i = 0; i < s.size();)
    {
        size_t skip = cli::skip_ansi_escape(s, i);
        if (skip != i)  // consumed an ANSI CSI escape, contributes no width
        {
            i = skip;
            continue;
        }
        if ((static_cast<unsigned char>(s[i]) & 0xc0U) != 0x80U)
        {  // count UTF-8 lead bytes, skip continuations
            ++width;
        }
        ++i;
    }
    return width;
}

auto pad(std::string_view cell, size_t width, Align align) -> std::string
{
    size_t visible = display_width(cell);
    size_t fill = width > visible ? width - visible : 0;
    if (align == Align::right)
    {
        return std::string(fill, ' ') + std::string(cell);
    }
    return std::string(cell) + std::string(fill, ' ');
}

auto rstrip(std::string s) -> std::string
{
    size_t end = s.find_last_not_of(' ');
    if (end == std::string::npos)
    {
        return "";
    }
    s.erase(end + 1);
    return s;
}

}  // namespace

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

auto Table::column(std::string_view header, Align align, size_t min_width)
    -> Table&
{
    columns_.push_back({std::string(header), align, min_width});
    return *this;
}

auto Table::row(std::vector<std::string> cells) -> Table&
{
    rows_.push_back({std::move(cells), false});
    return *this;
}

auto Table::rule() -> Table&
{
    rows_.push_back({{}, true});
    return *this;
}

auto Table::column_widths(bool include_rows) const -> std::vector<size_t>
{
    std::vector<size_t> widths(columns_.size());
    for (size_t i = 0; i < columns_.size(); ++i)
    {
        widths[i] = std::max(
            display_width(columns_[i].header), columns_[i].min_width);
    }
    if (include_rows)
    {
        for (const auto& r : rows_)
        {
            if (r.is_rule)
            {
                continue;
            }
            for (size_t i = 0; i < columns_.size() && i < r.cells.size(); ++i)
            {
                widths[i] = std::max(widths[i], display_width(r.cells[i]));
            }
        }
    }
    return widths;
}

auto Table::compose_row(
    std::span<const std::string> cells,
    std::span<const size_t> widths) const -> std::string
{
    std::string line(indent);
    for (size_t i = 0; i < widths.size(); ++i)
    {
        if (i != 0)
        {
            line += gap;
        }
        std::string_view cell = i < cells.size() ? std::string_view(cells[i])
                                                 : std::string_view();
        line += pad(cell, widths[i], columns_[i].align);
    }
    return rstrip(std::move(line));
}

auto Table::rule_line(std::span<const size_t> widths) -> std::string
{
    size_t inner = gap.size() * (widths.empty() ? 0 : widths.size() - 1);
    for (size_t w : widths)
    {
        inner += w;
    }
    return std::string(indent) + cli::separator(inner);
}

auto Table::header_block(std::span<const size_t> widths) const -> std::string
{
    std::vector<std::string> headers(columns_.size());
    for (size_t i = 0; i < columns_.size(); ++i)
    {
        headers[i] = columns_[i].header;
    }
    return compose_row(headers, widths) + '\n' + rule_line(widths);
}

auto Table::render() const -> std::string
{
    auto widths = column_widths(true);
    std::string out = header_block(widths);
    for (const auto& r : rows_)
    {
        out += '\n';
        out += r.is_rule ? rule_line(widths) : compose_row(r.cells, widths);
    }
    return out;
}

auto Table::stream_header() const -> std::string
{
    return header_block(column_widths(false));
}

auto Table::stream_row(std::span<const std::string> cells) const -> std::string
{
    return compose_row(cells, column_widths(false));
}

}  // namespace cli
