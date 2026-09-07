// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_COMPLETION_COMPLETION_GENERATOR_H_
#define APPS_CLI_COMPLETION_COMPLETION_GENERATOR_H_

#include <string>

namespace CLI
{
class App;
}  // namespace CLI

namespace cli
{
// Emit a completion script by introspecting the registered App tree: every
// visible subcommand, its options' names/description, enum choices recovered
// from the option type name ("<mode>:{A,D,AD}"), and file completion for
// path-like options.
auto generate_fish_completion(const CLI::App& root) -> std::string;

auto generate_bash_completion(const CLI::App& root) -> std::string;
}  // namespace cli

#endif  // APPS_CLI_COMPLETION_COMPLETION_GENERATOR_H_
