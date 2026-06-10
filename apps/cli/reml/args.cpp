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

    cmd.description(
        "Perform REML variance component estimation using average information "
        "algorithm");

    cmd.add_option(
           "-p,--pheno", "Phenotype file (TSV format: FID, IID, trait1, ...)")
        ->group("I/O")
        ->type_name("<PHENOTYPE>")
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
    cmd.add_option("--rand", "Random effects (TSV: FID, IID, effect1, ...)")
        ->group("I/O");
    cmd.add_option("-o,--out", "Output file prefix")
        ->group("I/O")
        ->type_name("<OUT>")
        ->default_val(std::string{"gelex"});

    cmd.add_option("--pheno-col", "Phenotype column index (0-based)")
        ->group("Model")
        ->default_val(2);

    cmd.add_option("--max-iter", "Max iterations")
        ->group("Runtime")
        ->default_val(100);
    cmd.add_option("--tol", "Convergence tolerance")
        ->group("Runtime")
        ->default_val(1e-6);
    cmd.add_option("-t,--threads", "CPU threads")
        ->group("Runtime")
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
