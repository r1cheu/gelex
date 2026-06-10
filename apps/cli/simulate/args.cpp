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

    cmd.formatter(cli::make_cli_formatter());
    cmd.description(
        "Simulate phenotypes based on genetic data and specified parameters");

    cmd.add_option("-b,--bfile", "PLINK binary file prefix (.bed/.bim/.fam)")
        ->group("I/O")
        ->type_name("<BFILE>")
        ->required();
    cmd.add_option("-o,--out", "Output file prefix")
        ->group("I/O")
        ->type_name("<OUT>")
        ->default_val(std::string{"sim.phen"});

    cmd.add_option(
           "--h2", "Narrow-sense heritability (0, 1); omit to disable additive")
        ->group("Model");
    cmd.add_option("--add-var", "Variances for additive effect classes")
        ->group("Model")
        ->type_name("<VARIANCES>")
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option(
           "--add-n",
           "SNP counts for additive effect classes (must match --add-var "
           "length)")
        ->group("Model")
        ->type_name("<COUNTS>")
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option(
           "--d2", "Dominance variance proportion (0, 1); h2+d2<1 in AD mode")
        ->group("Model");
    cmd.add_option("--dom-var", "Variances for dominance effect classes")
        ->group("Model")
        ->type_name("<VARIANCES>")
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option(
           "--dom-n",
           "SNP counts for dominance effect classes (must match --dom-var "
           "length)")
        ->group("Model")
        ->type_name("<COUNTS>")
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option(
           "--dom-pos-prob",
           "Probability of positive dominance effects; enables "
           "truncated-normal sampling [0, 1]")
        ->group("Model")
        ->type_name("<PROB>");
    cmd.add_option(
           "--geno-method,--gm",
           "Genotype processing: SH, CH, OSH, OCH, S, C, OS, OC, NS, NC "
           "(prefix: O=orth, N=NOIA; suffix: H=HWE-based)")
        ->group("Model")
        ->type_name("<STR>")
        ->default_val(std::string{"OS"});
    cmd.add_option("--seed", "Random seed")->group("Runtime")->default_val(42);

    cmd.footer(
        "Example:\n"
        "  gelex simulate -b geno\n\n"
        "Docs:\n"
        "  https://gelex.readthedocs.io/en/latest/cli/simulate.html");

    cmd.callback(
        [subcommand, &exit_code]()
        {
            exit_code = cli::execute_cli_command(*subcommand, simulate_execute);
        });
}
