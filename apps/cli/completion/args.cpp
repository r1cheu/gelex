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

#include "args.h"

#include <CLI/CLI.hpp>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "gelex/exception.h"

#include "cli/completion/completion_generator.h"

namespace
{
auto env_dir(const char* var, std::string_view fallback) -> std::string
{
    if (const char* value = std::getenv(var);
        value != nullptr && *value != '\0')
    {
        return value;
    }
    const char* home = std::getenv("HOME");
    if (home == nullptr || *home == '\0')
    {
        throw gelex::GelexException("cannot locate HOME for --install");
    }
    return std::string{home} + std::string{fallback};
}

auto run_completion(
    const CLI::App& root,
    const std::string& shell,
    bool install) -> void
{
    const std::string script = shell == "fish"
                                   ? cli::generate_fish_completion(root)
                                   : cli::generate_bash_completion(root);
    if (!install)
    {
        std::cout << script;
        return;
    }
    const auto& program = root.get_name();
    const std::filesystem::path path
        = shell == "fish"
              ? std::filesystem::path{env_dir("XDG_CONFIG_HOME", "/.config")}
                    / "fish" / "completions" / (program + ".fish")
              : std::filesystem::path{env_dir("XDG_DATA_HOME", "/.local/share")}
                    / "bash-completion" / "completions" / program;
    std::filesystem::create_directories(path.parent_path());
    std::ofstream out{path};
    if (!out)
    {
        throw gelex::GelexException(
            "cannot write completion script to " + path.string());
    }
    out << script;
    std::cerr << "Installed " << shell << " completion to " << path.string()
              << '\n';
}
}  // namespace

auto setup_completion_command(CLI::App& program, int& exit_code) -> void
{
    auto shell = std::make_shared<std::string>();
    auto install = std::make_shared<bool>(false);

    auto* subcommand = program.add_subcommand(
        "completion", "Generate a shell completion script");
    auto& cmd = *subcommand;
    cmd.group("");

    cmd.add_option("shell", *shell, "Target shell")
        ->required()
        ->check(CLI::IsMember(std::vector<std::string>{"fish", "bash"}));
    cmd.add_flag(
        "--install",
        *install,
        "Write into the shell's user completion directory instead of stdout");

    cmd.callback(
        [&program, &exit_code, shell, install]()
        {
            try
            {
                run_completion(program, *shell, *install);
                exit_code = 0;
            }
            catch (const std::exception& err)
            {
                std::cerr << "Error: " << err.what() << '\n';
                exit_code = 1;
            }
        });
}
