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

#include <CLI/CLI.hpp>
#include <fmt/format.h>
#include <iostream>
#include <string>

#include "cli/assoc/args.h"
#include "cli/color_formatter.h"
#include "cli/completion/args.h"
#include "cli/grm/args.h"
#include "cli/mcmc/args.h"
#include "cli/post/args.h"
#include "cli/predict/args.h"
#include "cli/reml/args.h"
#include "cli/simulate/args.h"
#include "version.h"

auto main(int argc, char* argv[]) -> int
{
    CLI::App program{
        fmt::format(
            R"(Gelex [version {}]
Genomic prediction, association testing, and variance-component estimation.

Issues and feature requests: https://github.com/r1cheu/gelex/issues
Docs: https://gelex.readthedocs.io/en/latest/)",
            PROJECT_VERSION),
        PROJECT_NAME};
    program.formatter(cli::make_cli_formatter());
    program.failure_message(cli::format_parse_error);
    program.set_version_flag("-v,--version", PROJECT_VERSION);
    program.require_subcommand(1);

    int exit_code = 0;
    setup_mcmc_command(program, exit_code);
    setup_simulate_command(program, exit_code);
    setup_predict_command(program, exit_code);
    setup_grm_command(program, exit_code);
    setup_assoc_command(program, exit_code);
    setup_post_command(program, exit_code);
    setup_reml_command(program, exit_code);
    setup_completion_command(program, exit_code);

    try
    {
        program.parse(argc, argv);
    }
    catch (const CLI::ParseError& err)
    {
        const CLI::App* active = &program;
        if (const auto entered = program.get_subcommands(); !entered.empty())
        {
            active = entered.back();
        }

        const auto& name = err.get_name();
        if (name == "CallForAllHelp")
        {
            std::cout << active->help("", CLI::AppFormatMode::All);
            return 0;
        }

        int command_tokens = 0;
        for (const auto* current = active; current != nullptr;
             current = current->get_parent())
        {
            command_tokens += 1;
        }
        // clap-style arg_required_else_help: a (sub)command invoked with no
        // arguments shows its full help instead of a missing-required error.
        if (name == "CallForHelp" || argc <= command_tokens)
        {
            std::cout << active->help();
            return 0;
        }

        return program.exit(err);
    }

    return exit_code;
}
