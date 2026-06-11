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
#include <vector>

#include <CLI/CLI.hpp>

#include "cli/cli_helper.h"
#include "command.h"

auto setup_assoc_command(CLI::App& program, int& exit_code) -> void
{
    auto* subcommand = program.add_subcommand("assoc");
    auto& cmd = *subcommand;

    cmd.description("Run mixed-model association testing");

    cli::add_common_io_options(cmd);
    cmd.add_option("-b,--bfile", "PLINK prefix; reads <prefix>.bed/.bim/.fam")
        ->group("I/O")
        ->type_name("<BFILE>")
        ->required();
    cmd.add_option(
           "--grm", "GRM prefix(es) without suffix; reads <prefix>.bin/.id")
        ->group("I/O")
        ->type_name("<GRM>")
        ->expected(1, -1)
        ->allow_extra_args()
        ->required();
    cmd.add_option("-o,--out", "Output prefix for .gwas.tsv")
        ->group("I/O")
        ->type_name("<OUT>")
        ->default_val(std::string{"gelex"});

    cmd.add_option("--transform", "Phenotype transform: none, dint, iint")
        ->group("Model")
        ->type_name("<TRANSFORM>")
        ->default_val(std::string{"none"})
        ->check(
            CLI::IsMember(std::vector<std::string>{"none", "dint", "iint"}));
    cmd.add_option("--int-offset", "Rank-INT offset k")
        ->group("Model")
        ->type_name("<K>")
        ->default_val(3.0 / 8.0);

    cmd.add_option("--test", "Wald test: single (df=1), joint (add+dom, df=2)")
        ->group("Model")
        ->type_name("<TEST>")
        ->default_val(std::string{"single"})
        ->check(CLI::IsMember(std::vector<std::string>{"single", "joint"}));
    cmd.add_option("--mode", "Effect mode for --test=single: A or D")
        ->group("Model")
        ->type_name("<MODE>")
        ->default_val(std::string{"A"})
        ->check(CLI::IsMember(std::vector<std::string>{"A", "D"}));
    cmd.add_flag("--loco", "Use leave-one-chromosome-out GRMs")->group("Model");
    cmd.add_option(
           "--geno-method,--gm",
           "Genotype coding: SH, CH, OSH, OCH, S, C, OS, OC, NS, NC")
        ->group("Model")
        ->type_name("<CODING>")
        ->default_val(std::string{"OCH"});

    cmd.add_option("--max-iter", "Maximum model-fit iterations")
        ->group("Runtime")
        ->type_name("<N>")
        ->default_val(100);
    cmd.add_option("--tol", "Convergence tolerance")
        ->group("Runtime")
        ->type_name("<TOL>")
        ->default_val(1e-6);
    cmd.add_option("-c,--chunk-size", "SNPs per chunk")
        ->group("Runtime")
        ->type_name("<N>")
        ->default_val(10000);
    cmd.add_option("-t,--threads", "CPU threads")
        ->group("Runtime")
        ->type_name("<N>")
        ->default_val(
            std::max(
                1, static_cast<int>(std::thread::hardware_concurrency() / 2)));

    cmd.footer(
        "Docs:\n"
        "  https://gelex.readthedocs.io/en/latest/cli/assoc.html");

    cmd.callback(
        [subcommand, &exit_code]()
        { exit_code = cli::execute_cli_command(*subcommand, assoc_execute); });
}
