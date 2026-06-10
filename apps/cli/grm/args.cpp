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

    cmd.description(
        "Compute genomic relationship matrix (GRM) from PLINK "
        "binary files and output in GCTA format");

    cmd.add_option("-b,--bfile", "PLINK binary file prefix (.bed/.bim/.fam)")
        ->group("I/O")
        ->type_name("<BFILE>")
        ->required();
    cmd.add_option("-o,--out", "Output file prefix")
        ->group("I/O")
        ->type_name("<OUT>")
        ->default_val(std::string{"grm"});

    cmd.add_option(
           "--geno-method,--gm",
           "Genotype coding: StandardizeHWE(SH), CenterHWE(CH),"
           " OrthStandardizeHWE(OSH), OrthCenterHWE(OCH),"
           " Standardize(S), Center(C), OrthStandardize(OS), OrthCenter(OC),"
           " NOIAStandardize(NS), NOIACenter(NC)")
        ->group("Model")
        ->type_name("<STR>")
        ->default_val(std::string{"OSH"});
    cmd.add_flag("--add", "Additive GRM")->group("Model");
    cmd.add_flag("--dom", "Dominance GRM")->group("Model");
    cmd.add_flag("--loco", "GRM per chromosome")->group("Model");

    cmd.add_option(
           "-c,--chunk-size", "SNPs per chunk for memory-efficient computation")
        ->group("Runtime")
        ->type_name("<SIZE>")
        ->default_val(10000);
    cmd.add_option("-t,--threads", "Threads (-1 for all cores)")
        ->group("Runtime")
        ->type_name("<N>")
        ->default_val(
            static_cast<int>(std::thread::hardware_concurrency() / 2));

    cmd.footer(
        "Example:\n"
        "  gelex grm -b geno --add\n\n"
        "Docs:\n"
        "  https://gelex.readthedocs.io/en/latest/cli/grm.html");

    cmd.callback(
        [subcommand, &exit_code]()
        { exit_code = cli::execute_cli_command(*subcommand, grm_execute); });
}
