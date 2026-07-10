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

#include <catch2/catch_test_macros.hpp>
#include <string>

#include "cli/table.h"

namespace
{

auto repeat(std::string_view unit, size_t n) -> std::string
{
    std::string out;
    for (size_t i = 0; i < n; ++i)
    {
        out += unit;
    }
    return out;
}

auto strip_ansi(std::string_view s) -> std::string
{
    std::string out;
    for (size_t i = 0; i < s.size();)
    {
        if (s[i] == 0x1b)
        {
            while (i < s.size() && s[i] != 'm')
            {
                ++i;
            }
            if (i < s.size())
            {
                ++i;
            }
            continue;
        }
        out += s[i++];
    }
    return out;
}

}  // namespace

TEST_CASE(
    "Table sizes columns to content and matches separator",
    "[cli][table]")
{
    cli::Table t;
    t.column("Chr", cli::Align::right);
    t.column("LogL", cli::Align::right);
    t.row({"1", "-10.5"});
    t.row({"2", "-9.25"});

    // Chr width = max(3, 1) = 3; LogL width = max(4, 5) = 5.
    std::string expected = "  Chr    LogL\n";
    expected += "  " + repeat("─", 11) + "\n";
    expected += "    1   -10.5\n";
    expected += "    2   -9.25";

    REQUIRE(t.render() == expected);
}

TEST_CASE("Table left-aligns and drops trailing padding", "[cli][table]")
{
    cli::Table t;
    t.column("Component", cli::Align::left);
    t.column("Estimate", cli::Align::right);
    t.row({"add", "44.446"});
    t.row({"Residual", "0.000"});

    // Component width = max(9, 3, 8) = 9; Estimate width = max(8, 6, 5) = 8.
    std::string expected = "  Component   Estimate\n";
    expected += "  " + repeat("─", 20) + "\n";
    expected += "  add" + repeat(" ", 11) + "44.446\n";
    expected += "  Residual" + repeat(" ", 7) + "0.000";

    REQUIRE(t.render() == expected);
}

TEST_CASE("Table rule() spans the full table width", "[cli][table]")
{
    cli::Table t;
    t.column("Chr", cli::Align::right);
    t.column("LogL", cli::Align::right);
    t.row({"1", "-10.5"});
    t.rule();
    t.row({"Mean", "-9.9"});

    // "Mean" widens the Chr column to 4; LogL stays 5.
    std::string rule = "  " + repeat("─", 12);
    std::string expected = "   Chr    LogL\n";
    expected += rule + "\n";
    expected += repeat(" ", 5) + "1   -10.5\n";
    expected += rule + "\n";
    expected += "  Mean    -9.9";

    REQUIRE(t.render() == expected);
}

TEST_CASE(
    "Table measures display width, not bytes, for styled cells",
    "[cli][table]")
{
    cli::Table t;
    t.column("V", cli::Align::right);
    // Rebecca-purple "AB": visible width 2 despite ANSI + escape bytes.
    t.row({"\x1b[38;2;102;51;153mAB\x1b[0m"});

    std::string out = t.render();
    auto newline = out.find('\n');
    std::string header = out.substr(0, newline);
    std::string data = out.substr(out.find_last_of('\n') + 1);

    // Column width = max("V"=1, "AB"=2) = 2, so both rows align to width 4.
    REQUIRE(strip_ansi(header).size() == 4);
    REQUIRE(strip_ansi(data).size() == 4);
}

TEST_CASE(
    "Table stream mode sizes columns from header and min_width",
    "[cli][table]")
{
    cli::Table t;
    t.column("Iter", cli::Align::right, 4);
    t.column("LogL", cli::Align::right, 10);
    t.column("V(2839.test.add)", cli::Align::right, 10);

    std::string header = t.stream_header();
    // Widths: 4, 10, 16 (label wider than min). Separator inner = 4+10+16+6=36.
    std::string expected_header = "  Iter         LogL   V(2839.test.add)\n";
    expected_header += "  " + repeat("─", 36);
    REQUIRE(header == expected_header);

    std::vector<std::string> cells{"1", "-8471.45", "32.03"};
    REQUIRE(t.stream_row(cells) == "     1     -8471.45              32.03");
}
