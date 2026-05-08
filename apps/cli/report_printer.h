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

#ifndef GELEX_CLI_REPORT_PRINTER_H_
#define GELEX_CLI_REPORT_PRINTER_H_

#include <string>
#include <string_view>
#include <utility>

#include <fmt/format.h>

namespace gelex::cli
{

// Owns the rule "every block-level output is preceded by exactly one blank
// line, every block sequence is gap-free unless punctuated by a block".
//
// OCP: only line/block/warn primitives are exposed; new section types are
// produced by free functions in formatter.h that return strings, then printed
// via block(...).
class ReportPrinter
{
   public:
    template <typename... Args>
    auto line(fmt::format_string<Args...> fmt_str, Args&&... args) -> void
    {
        emit_info(fmt::format(fmt_str, std::forward<Args>(args)...));
    }
    auto line(std::string_view msg) -> void { emit_info(std::string(msg)); }

    template <typename... Args>
    auto warn(fmt::format_string<Args...> fmt_str, Args&&... args) -> void
    {
        emit_warn(fmt::format(fmt_str, std::forward<Args>(args)...));
    }
    auto warn(std::string_view msg) -> void { emit_warn(std::string(msg)); }

    template <typename... Args>
    auto block(fmt::format_string<Args...> fmt_str, Args&&... args) -> void
    {
        ensure_blank();
        line(fmt_str, std::forward<Args>(args)...);
    }
    auto block(std::string_view msg) -> void
    {
        ensure_blank();
        line(msg);
    }

    auto ensure_blank() -> void;
    auto on_progress_finished() -> void { has_blank_ = false; }

   private:
    auto emit_info(std::string msg) -> void;
    auto emit_warn(std::string msg) -> void;

    bool has_blank_ = true;
};

auto printer() -> ReportPrinter&;

}  // namespace gelex::cli

#endif  // GELEX_CLI_REPORT_PRINTER_H_
