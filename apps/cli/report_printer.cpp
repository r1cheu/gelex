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

#include "report_printer.h"

#include <cstddef>
#include <string_view>

#include "cli/logging.h"

namespace cli
{

namespace
{

// One logger record must be one physical line: the file sink prefixes each
// record with "[time] [level]", so a multi-line payload would leave every line
// after the first without a prefix. Split on '\n' and emit each line
// separately.
template <typename Emit>
void emit_lines(std::string_view msg, Emit emit)
{
    std::size_t start = 0;
    while (true)
    {
        std::size_t nl = msg.find('\n', start);
        if (nl == std::string_view::npos)
        {
            emit(msg.substr(start));
            return;
        }
        emit(msg.substr(start, nl - start));
        start = nl + 1;
    }
}

}  // namespace

auto ReportPrinter::ensure_blank() -> void
{
    if (!has_blank_)
    {
        cli::logging::get()->info("");
        has_blank_ = true;
    }
}

auto ReportPrinter::emit_info(std::string_view msg) -> void
{
    emit_lines(
        msg,
        [](std::string_view line) { cli::logging::get()->info("{}", line); });
    has_blank_ = false;
}

auto ReportPrinter::emit_warn(std::string_view msg) -> void
{
    emit_lines(
        msg,
        [](std::string_view line) { cli::logging::get()->warn("{}", line); });
    has_blank_ = false;
}

auto printer() -> ReportPrinter&
{
    static ReportPrinter instance;
    return instance;
}

}  // namespace cli
