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

#include "color_formatter.h"

#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <CLI/CLI.hpp>

#include "cli_helper.h"

namespace cli
{

namespace
{

constexpr std::string_view RESET = "\033[0m";
constexpr std::string_view BOLD = "\033[1m";
constexpr std::string_view GREEN = "\033[32m";
constexpr std::string_view CYAN = "\033[36m";
constexpr std::string_view YELLOW = "\033[33m";

auto colored(std::string text, std::string_view color) -> std::string
{
    if (text.empty() || !is_tty())
    {
        return text;
    }
    return std::string{color} + text + std::string{RESET};
}

class ColorFormatter : public CLI::Formatter
{
   public:
    auto make_group(
        std::string group,
        bool is_positional,
        std::vector<const CLI::Option*> options) const -> std::string override
    {
        return CLI::Formatter::make_group(
            colored(std::move(group), std::string{BOLD} + std::string{GREEN}),
            is_positional,
            std::move(options));
    }

    auto make_option_name(const CLI::Option* option, bool is_positional) const
        -> std::string override
    {
        return colored(
            CLI::Formatter::make_option_name(option, is_positional), CYAN);
    }

    auto make_option_opts(const CLI::Option* option) const
        -> std::string override
    {
        return colored(CLI::Formatter::make_option_opts(option), YELLOW);
    }

    auto make_footer(const CLI::App* app) const -> std::string override
    {
        std::istringstream lines(app->get_footer());
        std::ostringstream out;
        std::string line;
        while (std::getline(lines, line))
        {
            if (line.ends_with(':'))
            {
                out << colored(
                    std::move(line), std::string{BOLD} + std::string{GREEN});
            }
            else
            {
                out << line;
            }
            out << '\n';
        }
        return out.str();
    }
};

}  // namespace

auto make_cli_formatter() -> std::shared_ptr<CLI::FormatterBase>
{
    auto formatter = std::make_shared<ColorFormatter>();
    formatter->enable_footer_formatting(false);
    return formatter;
}

}  // namespace cli
