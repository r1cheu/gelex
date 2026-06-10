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

#include "cli/cli_helper.h"
#include "command.h"

auto setup_predict_command(CLI::App& program, int& exit_code) -> void
{
    auto* subcommand = program.add_subcommand("predict");
    auto& cmd = *subcommand;

    cmd.description("Generate genomic predictions using fitted SNP effects");

    cmd.add_option("-b,--bfile", "PLINK binary prefix (.bed/.bim/.fam)")
        ->group("I/O")
        ->type_name("<BFILE>")
        ->required();
    cmd.add_option("-g,--gfile", "Fitted model prefix (.snpeff, .sbin, .param)")
        ->group("I/O")
        ->type_name("<GFILE>")
        ->required();
    cmd.add_option(
           "--qcovar", "Quantitative covariates (TSV: FID, IID, covar1, ...)")
        ->group("I/O")
        ->type_name("<QCOVAR>");
    cmd.add_option(
           "--dcovar", "Discrete covariates (TSV: FID, IID, factor1, ...)")
        ->group("I/O")
        ->type_name("<DCOVAR>");
    cmd.add_option("-o,--out", "Output file prefix")
        ->group("I/O")
        ->type_name("<OUT>")
        ->required();

    cmd.add_option("-c,--chunk-size", "SNPs per chunk")
        ->group("Runtime")
        ->default_val(10000);

    cmd.footer(
        "Example:\n"
        "  gelex predict -b geno -g model -o pred.tsv\n\n"
        "Docs:\n"
        "  https://gelex.readthedocs.io/en/latest/cli/predict.html");

    cmd.callback(
        [subcommand, &exit_code]()
        {
            exit_code = cli::execute_cli_command(*subcommand, predict_execute);
        });
}
