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

auto setup_mcmc_command(CLI::App& program, int& exit_code) -> void
{
    auto* subcommand = program.add_subcommand("mcmc");
    auto& cmd = *subcommand;

    cmd.description("Fit genomic prediction models with MCMC");

    cli::add_common_io_options(cmd);
    cmd.add_option("-b,--bfile", "PLINK prefix; reads <prefix>.bed/.bim/.fam")
        ->group("I/O")
        ->type_name("<BFILE>")
        ->required();
    cmd.add_option("-o,--out", "Output prefix for model files")
        ->group("I/O")
        ->type_name("<OUT>")
        ->default_val(std::string{"gelex"});

    cmd.add_option(
           "--geno-method,--gm",
           "Genotype coding: SH, CH, OSH, OCH, S, C, OS, OC, NS, NC")
        ->group("Model")
        ->type_name("<CODING>")
        ->default_val(std::string{"OSH"})
        ->check(cli::genotype_method_validator());
    cmd.add_option("-m,--method", "Bayesian model: RR, A, B, C, R, CD")
        ->group("Model")
        ->type_name("<MODEL>")
        ->default_val(std::string{"RR"})
        ->check(
            CLI::IsMember(
                std::vector<std::string>{"A", "B", "C", "R", "RR", "CD"}))
        ->required();
    cmd.add_option("--mode", "Effect mode: A, D, AD")
        ->group("Model")
        ->type_name("<MODE>")
        ->default_val(std::string{"A"})
        ->check(CLI::IsMember(std::vector<std::string>{"A", "D", "AD"}));
    cmd.add_option(
           "--random-pve", "Non-SNP random-effect variance fraction (0,1)")
        ->group("Model")
        ->type_name("<P>")
        ->check(cli::open_unit_interval());
    cmd.add_option("--h2", "Additive heritability (0,1)")
        ->group("Model")
        ->type_name("<P>")
        ->check(cli::open_unit_interval());
    cmd.add_option("--d2", "Dominance heritability (0,1)")
        ->group("Model")
        ->type_name("<P>")
        ->check(cli::open_unit_interval());
    cmd.add_option(
           "--dom-pos-prob",
           "Initial probability dominance effects are positive")
        ->group("Model")
        ->type_name("<P>")
        ->check(cli::open_unit_interval());
    cmd.add_option("--pi", "Additive mixture proportions (B/C/R)")
        ->group("Model")
        ->type_name("<P>")
        ->check(CLI::PositiveNumber)
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option("--dpi", "Dominance mixture proportions (B/C/R)")
        ->group("Model")
        ->type_name("<P>")
        ->check(CLI::PositiveNumber)
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option("--scale", "Additive variance multipliers (R)")
        ->group("Model")
        ->type_name("<SCALE>")
        ->check(cli::non_negative_number())
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option("--dscale", "Dominance variance multipliers (R)")
        ->group("Model")
        ->type_name("<SCALE>")
        ->check(cli::non_negative_number())
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option("--jpi", "Joint allocation proportions (CD)")
        ->group("Model")
        ->type_name("<P>")
        ->check(CLI::PositiveNumber)
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_flag("--sample-pi", "Sample additive mixture proportions")
        ->group("Model");
    cmd.add_flag("--sample-dpi", "Sample dominance mixture proportions")
        ->group("Model");
    cmd.add_flag("--sample-jpi", "Sample joint allocation proportions (CD)")
        ->group("Model");

    cmd.add_option("--iters", "MCMC iterations")
        ->group("Runtime")
        ->type_name("<N>")
        ->check(CLI::PositiveNumber)
        ->default_val(5000);
    cmd.add_option("--burn-in", "Burn-in iterations")
        ->group("Runtime")
        ->type_name("<N>")
        ->check(cli::non_negative_number())
        ->default_val(3000);
    cmd.add_option("--thin", "Sample every N iterations")
        ->group("Runtime")
        ->type_name("<N>")
        ->check(CLI::PositiveNumber)
        ->default_val(1);
    cmd.add_option("--seed", "Random seed")
        ->group("Runtime")
        ->type_name("<N>")
        ->default_val(42);
    cmd.add_option(
           "--checkpoint-step", "Checkpoint every N iterations; omit for final")
        ->group("Runtime")
        ->type_name("<N>")
        ->check(cli::non_negative_number());
    cmd.add_option("--from-ckpt", "Resume from checkpoint file")
        ->group("Runtime")
        ->type_name("<CKPT>")
        ->check(CLI::ExistingFile);
    cmd.add_option("-c,--chunk-size", "SNPs per chunk")
        ->group("Runtime")
        ->type_name("<N>")
        ->check(CLI::PositiveNumber)
        ->default_val(10000);
    cmd.add_option("-t,--threads", "CPU threads")
        ->group("Runtime")
        ->type_name("<N>")
        ->check(cli::non_negative_number())
        ->default_val(
            std::max(
                1, static_cast<int>(std::thread::hardware_concurrency() / 2)));
    cmd.add_flag("--mmap", "Store genotype chunks as mmap files")
        ->group("Runtime");

    cmd.footer(
        "Docs:\n"
        "  https://gelex.readthedocs.io/en/latest/cli/mcmc.html");

    cmd.callback(
        [subcommand, &exit_code]()
        { exit_code = cli::execute_cli_command(*subcommand, mcmc_execute); });
}
