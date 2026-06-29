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

auto setup_reml_command(CLI::App& program, int& exit_code) -> void
{
    auto config = std::make_shared<cli::RemlConfig>();
    auto* subcommand = program.add_subcommand("reml");
    auto& cmd = *subcommand;

    cmd.description("Estimate variance components with REML (AI algorithm)");

    cli::add_common_io_options(cmd, config->base_data);
    cmd.add_option(
           "--grm",
           config->grm_prefixes,
           "GRM prefix without suffix; reads <prefix>.bin/.id")
        ->group("I/O")
        ->type_name("<GRM>")
        ->expected(1, -1)
        ->allow_extra_args()
        ->required();
    cmd.add_option(
           "--drand",
           config->drand_path,
           "Discrete random-effect TSV (FID, IID, factor columns); each factor "
           "column becomes a variance component via one-hot ZZ^T")
        ->group("I/O")
        ->type_name("<DRAND>")
        ->check(CLI::ExistingFile);
    cmd.add_option(
           "--qrand",
           config->qrand_paths,
           "Quantitative random-effect matrix TSV (FID, IID, value columns); "
           "each file forms one linear-kernel component ZZ^T")
        ->group("I/O")
        ->type_name("<QRAND>")
        ->expected(1, -1)
        ->allow_extra_args()
        ->check(CLI::ExistingFile);
    cmd.add_option(
           "-o,--out",
           config->out_prefix,
           "Output prefix for .summary and .effects")
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

    cmd.add_option("--max-iter", config->max_iter, "Maximum AI-REML iterations")
        ->group("Runtime")
        ->type_name("<N>")
        ->check(CLI::PositiveNumber)
        ->capture_default_str();
    cmd.add_option(
           "--tol",
           config->tolerance,
           "Relative tolerance for variance-component convergence")
        ->group("Runtime")
        ->type_name("<TOL>")
        ->check(CLI::PositiveNumber)
        ->capture_default_str();
    cmd.add_option("-t,--threads", config->threads, "CPU threads")
        ->group("Runtime")
        ->type_name("<N>")
        ->check(cli::non_negative_number())
        ->capture_default_str();

    cmd.footer(
        "Docs:\n"
        "  https://gelex.readthedocs.io/en/latest/cli/reml.html");

    cmd.callback(
        [&cmd, &exit_code, config]()
        {
            exit_code = cli::execute_cli_command(
                cmd,
                "REML Variance Component Estimation",
                [config]() { return reml_execute(*config); });
        });
}
