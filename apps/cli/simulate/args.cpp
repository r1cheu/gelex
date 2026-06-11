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

auto setup_simulate_command(CLI::App& program, int& exit_code) -> void
{
    auto* subcommand = program.add_subcommand("simulate");
    auto& cmd = *subcommand;

    cmd.description("Simulate phenotypes from PLINK genotypes");

    cmd.add_option("-b,--bfile", "PLINK prefix; reads <prefix>.bed/.bim/.fam")
        ->group("I/O")
        ->type_name("<BFILE>")
        ->required();
    cmd.add_option("-o,--out", "Output prefix for .phen and .causal")
        ->group("I/O")
        ->type_name("<OUT>")
        ->default_val(std::string{"sim.phen"});

    cmd.add_option("--h2", "Additive heritability (0,1)")
        ->group("Model")
        ->type_name("<P>");
    cmd.add_option("--add-var", "Variances for additive effect classes")
        ->group("Model")
        ->type_name("<VAR>")
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option("--add-n", "SNP counts for additive effect classes")
        ->group("Model")
        ->type_name("<N>")
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option("--d2", "Dominance heritability (0,1)")
        ->group("Model")
        ->type_name("<P>");
    cmd.add_option("--dom-var", "Variances for dominance effect classes")
        ->group("Model")
        ->type_name("<VAR>")
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option("--dom-n", "SNP counts for dominance effect classes")
        ->group("Model")
        ->type_name("<N>")
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option(
           "--dom-pos-prob", "Probability dominance effects are positive")
        ->group("Model")
        ->type_name("<P>");
    cmd.add_option(
           "--geno-method,--gm",
           "Genotype coding: SH, CH, OSH, OCH, S, C, OS, OC, NS, NC")
        ->group("Model")
        ->type_name("<CODING>")
        ->default_val(std::string{"OS"});
    cmd.add_option("--seed", "Random seed")
        ->group("Runtime")
        ->type_name("<N>")
        ->default_val(42);

    cmd.footer(
        "Docs:\n"
        "  https://gelex.readthedocs.io/en/latest/cli/simulate.html");

    cmd.callback(
        [subcommand, &exit_code]()
        {
            exit_code = cli::execute_cli_command(*subcommand, simulate_execute);
        });
}
