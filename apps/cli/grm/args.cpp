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

#include <memory>
#include <string>
#include <vector>

#include <CLI/CLI.hpp>

#include "cli/cli_helper.h"
#include "command.h"
#include "config.h"

auto setup_grm_command(CLI::App& program, int& exit_code) -> void
{
    auto config = std::make_shared<cli::GrmConfig>();
    auto* subcommand = program.add_subcommand("grm");
    auto& cmd = *subcommand;

    cmd.description("Build genomic relationship matrices from PLINK data");

    cmd.add_option(
           "-b,--bfile",
           config->bfile,
           "PLINK prefix; reads <prefix>.bed/.bim/.fam")
        ->group("I/O")
        ->type_name("<BFILE>")
        ->required();
    cmd.add_option(
           "-o,--out",
           config->out,
           "Output prefix; writes <prefix>.<effect>.bin/.id")
        ->group("I/O")
        ->type_name("<OUT>")
        ->capture_default_str();

    cmd.add_option(
           "--geno-method,--gm",
           config->geno_method,
           "Genotype coding: SH, CH, OSH, OCH, S, C, OS, OC, NS, NC")
        ->group("Model")
        ->type_name("<CODING>")
        ->default_str("OSH")
        ->check(cli::genotype_method_validator());

    cmd.add_option("--mode", config->mode, "Effect mode: A, D, AD")
        ->group("Model")
        ->type_name("<MODE>")
        ->default_str("A")
        ->check(CLI::IsMember(std::vector<std::string>{"A", "D", "AD"}));
    cmd.add_flag("--loco", config->loco, "Write one GRM per chromosome")
        ->group("Model");

    cmd.add_option("-c,--chunk-size", config->chunk_size, "SNPs per chunk")
        ->group("Runtime")
        ->type_name("<N>")
        ->check(CLI::PositiveNumber)
        ->capture_default_str();
    cmd.add_option("-t,--threads", config->threads, "CPU threads")
        ->group("Runtime")
        ->type_name("<N>")
        ->check(cli::non_negative_number())
        ->capture_default_str();

    cmd.footer(
        "Docs:\n"
        "  https://gelex.readthedocs.io/en/latest/cli/grm.html");

    cmd.callback(
        [&cmd, &exit_code, config]()
        {
            exit_code = cli::execute_cli_command(
                cmd,
                "GRM Computation",
                [config]() { return grm_execute(*config); });
        });
}
