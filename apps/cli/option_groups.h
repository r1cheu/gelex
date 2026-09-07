// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_OPTION_GROUPS_H_
#define APPS_CLI_OPTION_GROUPS_H_

namespace CLI
{
class App;
}  // namespace CLI

namespace cli
{
struct BaseDataConfig;
struct RandomDesignDataConfig;
struct RemlDataConfig;

auto add_common_io_options(CLI::App& cmd, BaseDataConfig& config) -> void;

auto add_random_design_options(CLI::App& cmd, RandomDesignDataConfig& config)
    -> void;

auto add_random_effect_options(CLI::App& cmd, RemlDataConfig& config) -> void;
}  // namespace cli

#endif  // APPS_CLI_OPTION_GROUPS_H_
