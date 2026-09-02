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
#include <memory>

#include "cli/command_harness.h"
#include "cli/validators.h"
#include "command.h"
#include "config.h"

auto setup_post_command(CLI::App& program, int& exit_code) -> void
{
    auto config = std::make_shared<cli::PostConfig>();
    auto* subcommand = program.add_subcommand("post");
    auto& cmd = *subcommand;

    cmd.description("Summarize posterior diagnostics from MCMC samples");

    cmd.add_option(
           "--in", config->in, "MCMC output prefix(es); reads <prefix>.draws")
        ->group("I/O")
        ->type_name("<PREFIX>")
        ->expected(1, -1)
        ->allow_extra_args()
        ->required();
    cmd.add_option("-o,--out", config->out, "Output prefix for logs")
        ->group("I/O")
        ->type_name("<PREFIX>")
        ->capture_default_str();

    cmd.add_option("--hdpi", config->hdpi, "HPDI mass")
        ->group("Model")
        ->type_name("<P>")
        ->check(cli::open_unit_interval())
        ->capture_default_str();

    cmd.footer(
        "Docs:\n"
        "  https://gelex.readthedocs.io/en/latest/");

    cmd.callback(
        [&cmd, &exit_code, config]()
        {
            exit_code = cli::execute_cli_command(
                cmd,
                "MCMC Posterior Analysis",
                [config]() { return post_execute(*config); });
        });
}
