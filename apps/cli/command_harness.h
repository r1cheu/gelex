// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_COMMAND_HARNESS_H_
#define APPS_CLI_COMMAND_HARNESS_H_

#include <functional>
#include <string_view>

namespace CLI
{
class App;
}  // namespace CLI

namespace cli
{
auto execute_cli_command(
    const CLI::App& cmd,
    std::string_view banner_title,
    const std::function<int()>& execute_fn) -> int;
}  // namespace cli

#endif  // APPS_CLI_COMMAND_HARNESS_H_
