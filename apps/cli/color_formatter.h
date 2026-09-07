// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_COLOR_FORMATTER_H_
#define APPS_CLI_COLOR_FORMATTER_H_

#include <memory>
#include <string>

namespace CLI
{
class App;
class Error;
class FormatterBase;
}  // namespace CLI

namespace cli
{

auto make_cli_formatter() -> std::shared_ptr<CLI::FormatterBase>;

auto format_parse_error(const CLI::App* app, const CLI::Error& err)
    -> std::string;

}  // namespace cli

#endif  // APPS_CLI_COLOR_FORMATTER_H_
