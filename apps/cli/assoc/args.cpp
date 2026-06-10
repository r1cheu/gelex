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

    cmd.description(
        "Perform genome-wide association study using mixed linear model");

    cmd.add_option(
           "-p,--pheno", "Phenotype file (TSV format: FID, IID, trait1, ...)")
        ->group("I/O")
        ->type_name("<PHENOTYPE>")
        ->required();
    cmd.add_option(
           "--pheno-col",
           "0-based phenotype column after FID/IID; first trait=0")
        ->group("I/O")
        ->type_name("<COL>")
        ->default_val(0);
    cmd.add_option("-b,--bfile", "PLINK binary file prefix (.bed/.bim/.fam)")
        ->group("I/O")
        ->type_name("<BFILE>")
        ->required();
    cmd.add_option("--grm", "GRM file prefix(es). Can specify multiple GRMs.")
        ->group("I/O")
        ->type_name("<GRM>")
        ->expected(1, -1)
        ->allow_extra_args()
        ->required();
    cmd.add_option(
           "--qcovar", "Quantitative covariates (TSV: FID, IID, covar1, ...)")
        ->group("I/O");
    cmd.add_option(
           "--dcovar", "Discrete covariates (TSV: FID, IID, factor1, ...)")
        ->group("I/O");
    cmd.add_option("-o,--out", "Output file prefix")
        ->group("I/O")
        ->type_name("<OUT>")
        ->default_val(std::string{"gelex"});

    cmd.add_option(
           "--transform",
           "Phenotype transformation: none, dint (Direct INT), iint (Indirect "
           "INT)")
        ->group("Model")
        ->type_name("<TRANSFORM>")
        ->default_val(std::string{"none"})
        ->check(
            CLI::IsMember(std::vector<std::string>{"none", "dint", "iint"}));
    cmd.add_option("--int-offset", "INT offset parameter k (Blom offset)")
        ->group("Model")
        ->default_val(3.0 / 8.0);

    cmd.add_option(
           "--test",
           "Wald test mode: single (one-effect, df=1), joint (add+dom, df=2)")
        ->group("Model")
        ->type_name("<TEST>")
        ->default_val(std::string{"single"})
        ->check(CLI::IsMember(std::vector<std::string>{"single", "joint"}));
    cmd.add_option(
           "--mode",
           "Association mode: A (additive), D (dominance). "
           "Ignored when --test=joint (always AD)")
        ->group("Model")
        ->type_name("<MODE>")
        ->default_val(std::string{"A"})
        ->check(CLI::IsMember(std::vector<std::string>{"A", "D"}));
    cmd.add_flag("--loco", "Leave-One-Chromosome-Out analysis")->group("Model");
    cmd.add_option(
           "--geno-method,--gm",
           "Genotype processing: SH, CH, OSH, OCH, S, C, OS, OC, NS, NC "
           "(prefix: O=orth, N=NOIA; suffix: H=HWE-based)")
        ->group("Model")
        ->type_name("<STR>")
        ->default_val(std::string{"OCH"});

    cmd.add_option("--max-iter", "Max iterations")
        ->group("Runtime")
        ->default_val(100);
    cmd.add_option("--tol", "Convergence tolerance")
        ->group("Runtime")
        ->default_val(1e-6);
    cmd.add_option("-c,--chunk-size", "SNPs per chunk")
        ->group("Runtime")
        ->default_val(10000);
    cmd.add_option("-t,--threads", "CPU threads")
        ->group("Runtime")
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
