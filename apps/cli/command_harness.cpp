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

#include "cli/command_harness.h"

#include <CLI/CLI.hpp>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fmt/base.h>
#include <fmt/format.h>
#include <functional>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include "cli/formatter.h"
#include "cli/logging.h"
#include "report_printer.h"
#include "theme.h"
#include "version.h"

namespace cli
{
namespace
{

auto report_command_line(const CLI::App& cmd) -> void
{
    auto& p = cli::printer();

    auto cwd = std::filesystem::current_path().string();
    if (const char* home = std::getenv("HOME");
        home != nullptr && cwd.starts_with(home))
    {
        cwd.replace(0, std::string_view{home}.size(), "~");
    }
    p.block(cli::section("Command line ({})", cwd));

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
            p.line("  {}", name);
        }
        else
        {
            constexpr size_t wrap_width = 80;
            constexpr std::string_view cont_indent = "    ";
            std::string line = fmt::format("  {}", name);
            bool has_value = false;
            for (const auto& value : values)
            {
                if (has_value && line.size() + 1 + value.size() > wrap_width)
                {
                    p.line("{}", line);
                    line = fmt::format("{}{}", cont_indent, value);
                }
                else
                {
                    line += ' ';
                    line += value;
                }
                has_value = true;
            }
            p.line("{}", line);
        }
    }
}

}  // namespace

auto execute_cli_command(
    const CLI::App& cmd,
    std::string_view banner_title,
    const std::function<int()>& execute_fn) -> int
{
    try
    {
        cli::logging::initialize(cmd.get_option("--out")->as<std::string>());
        cli::printer().block(
            cli::command_banner(PROJECT_VERSION, banner_title));
        report_command_line(cmd);
        auto start = std::chrono::steady_clock::now();
        auto result = execute_fn();
        auto elapsed = std::chrono::duration<double>(
                           std::chrono::steady_clock::now() - start)
                           .count();
        if (cli::logging::get())
        {
            cli::printer().block(cli::done_message(elapsed));
        }
        return result;
    }
    catch (const std::exception& e)
    {
        auto logger = cli::logging::get();
        if (logger)
        {
            logger->error("{}", e.what());
        }
        else
        {
            std::cerr << error_marker() << e.what() << "\n";
        }
        return 1;
    }
    catch (...)
    {
        auto logger = cli::logging::get();
        if (logger)
        {
            logger->error("unknown exception");
        }
        else
        {
            std::cerr << error_marker() << "unknown exception\n";
        }
        return 1;
    }
}

}  // namespace cli
