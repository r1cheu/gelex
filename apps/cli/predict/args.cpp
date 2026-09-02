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
#include "command.h"
#include "config.h"

auto setup_predict_command(CLI::App& program, int& exit_code) -> void
{
    auto config = std::make_shared<cli::PredictConfig>();
    auto* subcommand = program.add_subcommand("predict");
    auto& cmd = *subcommand;

    cmd.description("Predict phenotypes from fitted SNP effects");

    cmd.add_option(
           "-b,--bfile",
           config->bfile,
           "PLINK prefix; reads <prefix>.bed/.bim/.fam")
        ->group("I/O")
        ->type_name("<BFILE>")
        ->required();
    cmd.add_option(
           "-g,--gfile",
           config->gfile,
           "Fitted model prefix; reads <prefix>.snpeff/.snplut/.param")
        ->group("I/O")
        ->type_name("<GFILE>")
        ->required();
    cmd.add_option(
           "--qcovar",
           config->qcovar,
           "Quantitative covariate TSV with FID, IID, numeric columns")
        ->group("I/O")
        ->type_name("<QCOVAR>")
        ->check(CLI::ExistingFile);
    cmd.add_option(
           "--dcovar",
           config->dcovar,
           "Discrete covariate TSV with FID, IID, factor columns")
        ->group("I/O")
        ->type_name("<DCOVAR>")
        ->check(CLI::ExistingFile);
    cmd.add_option("-o,--out", config->out, "Output prefix for .pred.tsv")
        ->group("I/O")
        ->type_name("<OUT>")
        ->required();

    cmd.add_option("-c,--chunk-size", config->chunk_size, "SNPs per chunk")
        ->group("Runtime")
        ->type_name("<N>")
        ->check(CLI::PositiveNumber)
        ->capture_default_str();

    cmd.footer(
        "Docs:\n"
        "  https://gelex.readthedocs.io/en/latest/cli/predict.html");

    cmd.callback(
        [&cmd, &exit_code, config]()
        {
            exit_code = cli::execute_cli_command(
                cmd,
                "Genomic Prediction",
                [config]() { return predict_execute(*config); });
        });
}
