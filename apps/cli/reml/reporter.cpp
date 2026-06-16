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

#include "reporter.h"

#include <algorithm>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <CLI/CLI.hpp>

#include "cli/report_printer.h"
#include "gelex/infra/logging/formatter.h"
#include "version.h"

namespace cli
{

auto RemlCommandReporter::show_banner() const -> void
{
    cli::printer().block(
        gelex::command_banner(
            PROJECT_VERSION, "REML Variance Component Estimation"));
}

auto RemlCommandReporter::show_config(const CLI::App& cmd) const -> void
{
    auto& p = cli::printer();
    p.block(gelex::section("Options in effect:"));

    std::vector<const CLI::Option*> printed;
    for (const auto* option : cmd.parse_order())
    {
        if (option->count() == 0)
        {
            continue;
        }
        if (std::ranges::find(printed, option) != printed.end())
        {
            continue;
        }
        printed.push_back(option);

        std::string name;
        const auto& long_names = option->get_lnames();
        const auto& short_names = option->get_snames();
        if (!long_names.empty())
        {
            name = "--" + long_names.front();
        }
        else if (!short_names.empty())
        {
            name = "-" + short_names.front();
        }
        if (name.empty())
        {
            continue;
        }

        std::vector<std::string> values;
        for (const auto& value : option->results())
        {
            if (!value.empty() && value != "true")
            {
                values.push_back(value);
            }
        }
        if (values.empty())
        {
            p.line("  {}", gelex::rebecca_purple(name));
        }
        else
        {
            p.line(
                "  {} {}", gelex::rebecca_purple(name), fmt::join(values, " "));
        }
    }
}

}  // namespace cli
