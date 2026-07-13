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
// from the option type name ("<MODE>:{A,D,AD}"), and file completion for
// path-like options.
auto generate_fish_completion(const CLI::App& root) -> std::string;

auto generate_bash_completion(const CLI::App& root) -> std::string;
}  // namespace cli

#endif  // APPS_CLI_COMPLETION_COMPLETION_GENERATOR_H_
