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
#include <thread>

#include <CLI/CLI.hpp>

#include "cli/cli_helper.h"
#include "command.h"

auto setup_grm_command(CLI::App& program, int& exit_code) -> void
{
    auto* subcommand = program.add_subcommand("grm");
    auto& cmd = *subcommand;

    cmd.description("Build genomic relationship matrices from PLINK data");

    cmd.add_option("-b,--bfile", "PLINK prefix; reads <prefix>.bed/.bim/.fam")
        ->group("I/O")
        ->type_name("<BFILE>")
        ->required();
    cmd.add_option(
           "-o,--out", "Output prefix; writes <prefix>.<effect>.bin/.id")
        ->group("I/O")
        ->type_name("<OUT>")
        ->default_val(std::string{"grm"});

    cmd.add_option(
           "--geno-method,--gm",
           "Genotype coding: SH, CH, OSH, OCH, S, C, OS, OC, NS, NC")
        ->group("Model")
        ->type_name("<CODING>")
        ->default_val(std::string{"OSH"})
        ->check(cli::genotype_method_validator());
    cmd.add_flag("--add", "Write additive GRM")->group("Model");
    cmd.add_flag("--dom", "Write dominance GRM")->group("Model");
    cmd.add_flag("--loco", "Write one GRM per chromosome")->group("Model");

    cmd.add_option("-c,--chunk-size", "SNPs per chunk")
        ->group("Runtime")
        ->type_name("<N>")
        ->check(CLI::PositiveNumber)
        ->default_val(10000);
    cmd.add_option("-t,--threads", "CPU threads")
        ->group("Runtime")
        ->type_name("<N>")
        ->check(CLI::NonNegativeNumber)
        ->default_val(
            static_cast<int>(std::thread::hardware_concurrency() / 2));

    cmd.footer(
        "Docs:\n"
        "  https://gelex.readthedocs.io/en/latest/cli/grm.html");

    cmd.callback(
        [subcommand, &exit_code]()
        { exit_code = cli::execute_cli_command(*subcommand, grm_execute); });
}
