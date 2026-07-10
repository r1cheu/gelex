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

#ifndef APPS_CLI_TABLE_H_
#define APPS_CLI_TABLE_H_

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace cli
{

enum class Align : std::uint8_t
{
    left,
    right
};

std::string separator(size_t width = 70, const std::string& c = "─");
std::string table_separator(size_t width = 70);
std::string
named_section(std::string_view name, size_t width = 32, size_t indent = 0);

// Fixed-gap column layout that sizes columns to their visible content and draws
// separators matching the real table width. Cells may carry ANSI styling and
// multi-byte UTF-8: width is measured by display columns, not bytes, so styled
// cells stay aligned (which fixed-width `fmt` padding cannot do).
//
// Two modes share the same layout core:
//   static  — collect every row via row()/rule(), then render(); columns size
//   to
//             the widest of header, cells, and min_width.
//   stream  — call stream_header() once, then stream_row() per row as data
//             arrives; columns size to max(header, min_width) up front, so give
//             each column a min_width wide enough for its values.
class Table
{
   public:
    auto column(
        std::string_view header,
        Align align = Align::right,
        size_t min_width = 0) -> Table&;

    auto row(std::vector<std::string> cells) -> Table&;
    auto rule() -> Table&;
    auto render() const -> std::string;

    auto stream_header() const -> std::string;
    auto stream_row(std::span<const std::string> cells) const -> std::string;

   private:
    struct Column
    {
        std::string header;
        Align align;
        size_t min_width;
    };
    struct Row
    {
        std::vector<std::string> cells;
        bool is_rule;
    };

    auto column_widths(bool include_rows) const -> std::vector<size_t>;
    auto compose_row(
        std::span<const std::string> cells,
        std::span<const size_t> widths) const -> std::string;
    auto header_block(std::span<const size_t> widths) const -> std::string;
    static auto rule_line(std::span<const size_t> widths) -> std::string;

    std::vector<Column> columns_;
    std::vector<Row> rows_;
};

}  // namespace cli

#endif  // APPS_CLI_TABLE_H_
