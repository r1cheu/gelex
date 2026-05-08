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

#include <utility>

#include <spdlog/logger.h>

#include "gelex/infra/logger.h"

namespace gelex::cli
{

auto ReportPrinter::ensure_blank() -> void
{
    if (!has_blank_)
    {
        gelex::logging::get()->info("");
        has_blank_ = true;
    }
}

auto ReportPrinter::emit_info(std::string msg) -> void
{
    gelex::logging::get()->info("{}", std::move(msg));
    has_blank_ = false;
}

auto ReportPrinter::emit_warn(std::string msg) -> void
{
    gelex::logging::get()->warn("{}", std::move(msg));
    has_blank_ = false;
}

auto printer() -> ReportPrinter&
{
    static ReportPrinter instance;
    return instance;
}

}  // namespace gelex::cli
