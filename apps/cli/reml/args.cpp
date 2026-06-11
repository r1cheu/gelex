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

#include <algorithm>
#include <string>
#include <thread>

#include <CLI/CLI.hpp>

#include "cli/cli_helper.h"
#include "command.h"

auto setup_reml_command(CLI::App& program, int& exit_code) -> void
{
    auto* subcommand = program.add_subcommand("reml");
    auto& cmd = *subcommand;

    cmd.description("Estimate variance components with REML (AI algorithm)");

    cli::add_common_io_options(cmd);
    cmd.add_option("--grm", "GRM prefix without suffix; reads <prefix>.bin/.id")
        ->group("I/O")
        ->type_name("<GRM>")
        ->expected(1, -1)
        ->allow_extra_args()
        ->required();
    cmd.add_option(
           "--rand", "Random-effect factor TSV with FID, IID, factor columns")
        ->group("I/O")
        ->type_name("<RAND>");
    cmd.add_option("-o,--out", "Output prefix for .summary and .effects")
        ->group("I/O")
        ->type_name("<OUT>")
        ->default_val(std::string{"gelex"});

    cmd.add_option("--max-iter", "Maximum AI-REML iterations")
        ->group("Runtime")
        ->type_name("<N>")
        ->default_val(100);
    cmd.add_option(
           "--tol", "Relative tolerance for variance-component convergence")
        ->group("Runtime")
        ->type_name("<TOL>")
        ->default_val(1e-6);
    cmd.add_option("-t,--threads", "CPU threads")
        ->group("Runtime")
        ->type_name("<N>")
        ->default_val(
            std::max(
                1, static_cast<int>(std::thread::hardware_concurrency() / 2)));

    cmd.footer(
        "Docs:\n"
        "  https://gelex.readthedocs.io/en/latest/cli/reml.html");

    cmd.callback(
        [subcommand, &exit_code]()
        { exit_code = cli::execute_cli_command(*subcommand, reml_execute); });
}
