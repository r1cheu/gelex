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

auto setup_assoc_command(CLI::App& program, int& exit_code) -> void
{
    auto config = std::make_shared<cli::AssocConfig>();
    auto* subcommand = program.add_subcommand("assoc");
    auto& cmd = *subcommand;

    cmd.description("Run mixed-model association testing");

    cli::add_common_io_options(cmd, config->base_data);
    cmd.add_option(
           "-b,--bfile",
           config->bfile,
           "PLINK prefix; reads <prefix>.bed/.bim/.fam")
        ->group("I/O")
        ->type_name("<BFILE>")
        ->required();
    cmd.add_option(
           "--grm",
           config->grm,
           "GRM prefix(es) without suffix; reads <prefix>.bin/.id")
        ->group("I/O")
        ->type_name("<GRM>")
        ->expected(1, -1)
        ->allow_extra_args()
        ->required();
    cmd.add_option("-o,--out", config->out, "Output prefix for .gwas.tsv")
        ->group("I/O")
        ->type_name("<OUT>")
        ->capture_default_str();

    cmd.add_option(
           "--transform",
           config->base_data.transform,
           "Phenotype transform: none, dint, iint")
        ->group("Model")
        ->type_name("<TRANSFORM>")
        ->capture_default_str()
        ->check(
            CLI::IsMember(std::vector<std::string>{"none", "dint", "iint"}));
    cmd.add_option(
           "--int-offset", config->base_data.int_offset, "Rank-INT offset k")
        ->group("Model")
        ->type_name("<K>")
        ->capture_default_str();

    cmd.add_option(
           "--mode",
           config->mode,
           "Wald test: A, D (single, df=1); AD (joint add+dom, df=2)")
        ->group("Model")
        ->type_name("<MODE>")
        ->default_str("A")
        ->check(CLI::IsMember(std::vector<std::string>{"A", "D", "AD"}));
    cmd.add_flag("--loco", config->loco, "Use leave-one-chromosome-out GRMs")
        ->group("Model");
    cmd.add_option(
           "--geno-method,--gm",
           config->geno_method,
           "Genotype coding: SH, CH, OSH, OCH, S, C, OS, OC, NS, NC")
        ->group("Model")
        ->type_name("<CODING>")
        ->default_str("OCH")
        ->check(cli::genotype_method_validator());

    cmd.add_option(
           "--max-iter", config->max_iter, "Maximum model-fit iterations")
        ->group("Runtime")
        ->type_name("<N>")
        ->check(CLI::PositiveNumber)
        ->capture_default_str();
    cmd.add_option("--tol", config->tolerance, "Convergence tolerance")
        ->group("Runtime")
        ->type_name("<TOL>")
        ->check(CLI::PositiveNumber)
        ->capture_default_str();
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
        "  https://gelex.readthedocs.io/en/latest/cli/assoc.html");

    cmd.callback(
        [&cmd, &exit_code, config]()
        {
            exit_code = cli::execute_cli_command(
                cmd,
                "GWAS Analysis",
                [config]() { return assoc_execute(*config); });
        });
}
