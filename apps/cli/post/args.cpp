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

#include <string>

#include <CLI/CLI.hpp>

#include "cli/cli_helper.h"
#include "command.h"

auto setup_post_command(CLI::App& program, int& exit_code) -> void
{
    auto* subcommand = program.add_subcommand("post");
    auto& cmd = *subcommand;

    cmd.description("Summarize posterior diagnostics from MCMC samples");

    cmd.add_option("--in", "MCMC output prefix(es); reads <prefix>.samples")
        ->group("I/O")
        ->type_name("<PREFIX>")
        ->expected(1, -1)
        ->allow_extra_args()
        ->required();
    cmd.add_option("-g,--gfile", "Genotype binary prefix for SNP diagnostics")
        ->group("I/O")
        ->type_name("<GFILE>");
    cmd.add_option("-o,--out", "Output prefix for logs")
        ->group("I/O")
        ->type_name("<PREFIX>")
        ->default_val(std::string{"gelex_post"});

    cmd.add_option("--hdpi", "HPDI mass")
        ->group("Model")
        ->type_name("<P>")
        ->default_val(0.94);

    cmd.footer(
        "Docs:\n"
        "  https://gelex.readthedocs.io/en/latest/");

    cmd.callback(
        [subcommand, &exit_code]()
        { exit_code = cli::execute_cli_command(*subcommand, post_execute); });
}
