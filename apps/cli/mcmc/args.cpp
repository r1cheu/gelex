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
#include <string>

#include "gelex/genetic_mode.h"

#include "cli/command_harness.h"
#include "cli/lexical_cast.h"
#include "cli/option_groups.h"
#include "cli/validators.h"
#include "command.h"
#include "config.h"

auto setup_mcmc_command(CLI::App& program, int& exit_code) -> void
{
    auto config = std::make_shared<cli::McmcConfig>();
    auto* subcommand = program.add_subcommand("mcmc");
    auto& cmd = *subcommand;

    cmd.description("Fit genomic prediction models with MCMC");

    cli::add_common_io_options(cmd, config->base_data);
    cli::add_random_design_options(cmd, config->random);
    cmd.add_option(
           "-b,--bfile",
           config->bfile,
           "PLINK prefix; reads <prefix>.bed/.bim/.fam")
        ->group("I/O")
        ->type_name("<BFILE>")
        ->required();
    cmd.add_option("-o,--out", config->out, "Output prefix for model files")
        ->group("I/O")
        ->type_name("<OUT>")
        ->capture_default_str();

    cmd.add_option(
           "--geno-method,--gm",
           config->geno_method,
           "Genotype coding: SH, CH, OSH, OCH, S, C, OS, OC, NS, NC")
        ->group("Model")
        ->type_name("<CODING>")
        ->default_str("OSH")
        ->check(cli::genotype_method_validator());
    cmd.add_option_function<std::string>(
           "-m,--method",
           cli::lexical_assigner(config->method),
           "Bayesian model: RR, A, B, C, R, CD")
        ->group("Model")
        ->type_name("<MODEL>")
        ->default_str("RR")
        ->check(cli::bayes_method_validator());
    cmd.add_option_function<std::string>(
           "--mode",
           cli::lexical_assigner(config->mode),
           "Effect mode: A, D, AD")
        ->group("Model")
        ->type_name("<MODE>")
        ->default_str("A")
        ->check(cli::genetic_mode_set_validator());
    cmd.add_option(
           "--h2",
           config->genetic_variance_shares.get<gelex::GeneticMode::A>(),
           "Additive heritability (0,1)")
        ->group("Model")
        ->type_name("<P>")
        ->check(cli::open_unit_interval());
    cmd.add_option(
           "--d2",
           config->genetic_variance_shares.get<gelex::GeneticMode::D>(),
           "Dominance heritability (0,1)")
        ->group("Model")
        ->type_name("<P>")
        ->check(cli::open_unit_interval());
    cmd.add_option(
           "--dom-pos-prob",
           config->dominance_positive_probability,
           "Initial probability dominance effects are positive")
        ->group("Model")
        ->type_name("<P>")
        ->check(cli::open_unit_interval());
    cmd.add_option(
           "--random-pve",
           config->random_pve,
           "Total variance fraction for non-SNP random effects")
        ->group("Model")
        ->type_name("<P>")
        ->check(cli::open_unit_interval());
    cmd.add_option(
           "--pi",
           config->mixture_probabilities.get<gelex::GeneticMode::A>(),
           "Additive mixture proportions (B/C/R)")
        ->group("Model")
        ->type_name("<P>")
        ->check(CLI::PositiveNumber)
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option(
           "--dpi",
           config->mixture_probabilities.get<gelex::GeneticMode::D>(),
           "Dominance mixture proportions (B/C/R)")
        ->group("Model")
        ->type_name("<P>")
        ->check(CLI::PositiveNumber)
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option(
           "--scale",
           config->mixture_scales.get<gelex::GeneticMode::A>(),
           "Additive variance multipliers (R)")
        ->group("Model")
        ->type_name("<SCALE>")
        ->check(cli::non_negative_number())
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option(
           "--dscale",
           config->mixture_scales.get<gelex::GeneticMode::D>(),
           "Dominance variance multipliers (R)")
        ->group("Model")
        ->type_name("<SCALE>")
        ->check(cli::non_negative_number())
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option(
           "--jpi",
           config->mixture_probabilities.joint(),
           "Joint allocation proportions (CD)")
        ->group("Model")
        ->type_name("<P>")
        ->check(CLI::PositiveNumber)
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option("--iters", config->iters, "MCMC iterations")
        ->group("Runtime")
        ->type_name("<N>")
        ->check(CLI::PositiveNumber)
        ->capture_default_str();
    cmd.add_option("--burn-in", config->burn_in, "Burn-in iterations")
        ->group("Runtime")
        ->type_name("<N>")
        ->check(cli::non_negative_number())
        ->capture_default_str();
    cmd.add_option("--thin", config->thin, "Sample every N iterations")
        ->group("Runtime")
        ->type_name("<N>")
        ->check(CLI::PositiveNumber)
        ->capture_default_str();
    cmd.add_option("--seed", config->seed, "Random seed")
        ->group("Runtime")
        ->type_name("<N>")
        ->capture_default_str();
    cmd.add_option("-t,--threads", config->threads, "CPU threads")
        ->group("Runtime")
        ->type_name("<N>")
        ->check(cli::non_negative_number())
        ->capture_default_str();
    cmd.footer(
        "Docs:\n"
        "  https://gelex.readthedocs.io/en/latest/cli/mcmc.html");

    cmd.callback(
        [&cmd, &exit_code, config]()
        {
            exit_code = cli::execute_cli_command(
                cmd,
                "Model Fitting (MCMC)",
                [config]() { return mcmc_execute(*config); });
        });
}
