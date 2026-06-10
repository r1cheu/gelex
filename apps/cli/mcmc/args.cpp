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

    cmd.description(
        "Fit genomic prediction models using MCMC (Gibbs sampling)");

    cmd.add_option(
           "-p,--pheno", "Phenotype file (TSV format: FID, IID, trait1, ...)")
        ->group("I/O")
        ->type_name("<PHENOTYPE>")
        ->required();
    auto* bfile = cmd.add_option(
        "-b,--bfile", "PLINK binary file prefix (.bed/.bim/.fam)");
    bfile->group("I/O")->type_name("<BFILE>");
    auto* gfile
        = cmd.add_option("-g,--gfile", "Gelex binary genotype file prefix");
    gfile->group("I/O")->type_name("<GFILE>");
    bfile->excludes(gfile);
    gfile->excludes(bfile);
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

    cmd.add_option("--pheno-col", "Phenotype column index (0-based)")
        ->group("Model")
        ->default_val(2);
    cmd.add_option(
           "--geno-method,--gm",
           "Genotype method: StandardizeHWE(SH), CenterHWE(CH),"
           " OrthStandardizeHWE(OSH), OrthCenterHWE(OCH),"
           " Standardize(S), Center(C), OrthStandardize(OS), OrthCenter(OC),"
           " NOIAStandardize(NS), NOIACenter(NC)")
        ->group("Model")
        ->type_name("<STR>")
        ->default_val(std::string{"OSH"});
    cmd.add_option("-m,--method", "Bayesian method: RR, A, B, C, R, CD")
        ->group("Model")
        ->type_name("<METHOD>")
        ->default_val(std::string{"RR"})
        ->check(
            CLI::IsMember(
                std::vector<std::string>{"A", "B", "C", "R", "RR", "CD"}))
        ->required();
    cmd.add_option(
           "--mode",
           "Genetic effect mode: A (additive only), D (dominance only), "
           "AD (additive + dominance)")
        ->group("Model")
        ->type_name("<MODE>")
        ->default_val(std::string{"A"})
        ->check(CLI::IsMember(std::vector<std::string>{"A", "D", "AD"}));
    cmd.add_option(
           "--random-pve",
           "Variance fraction for non-SNP random effects (0, 1)")
        ->group("Model");
    cmd.add_option("--h2", "Additive heritability (0, 1)")->group("Model");
    cmd.add_option("--d2", "Dominance heritability (0, 1)")->group("Model");
    cmd.add_option(
           "--dom-pos-prob",
           "Initial probability that an active dominance effect is positive")
        ->group("Model");
    cmd.add_option("--pi", "Additive mixture proportions (B/C/R)")
        ->group("Model")
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option("--dpi", "Dominance mixture proportions (B/C/R)")
        ->group("Model")
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option("--scale", "Additive variance multipliers (R)")
        ->group("Model")
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option("--dscale", "Dominance variance multipliers (R)")
        ->group("Model")
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option("--jpi", "Joint allocation proportions (CD)")
        ->group("Model")
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_flag("--sample-pi", "Sample additive mixture proportions")
        ->group("Model");
    cmd.add_flag("--sample-dpi", "Sample dominance mixture proportions")
        ->group("Model");
    cmd.add_flag("--sample-jpi", "Sample joint allocation proportions (CD)")
        ->group("Model");

    cmd.add_option("--iters", "MCMC iterations to run")
        ->group("Runtime")
        ->default_val(5000);
    cmd.add_option("--burn-in", "Burn-in iterations to discard")
        ->group("Runtime")
        ->default_val(3000);
    cmd.add_option("--thin", "Thinning interval for samples")
        ->group("Runtime")
        ->default_val(1);
    cmd.add_option("--seed", "Random seed for MCMC")
        ->group("Runtime")
        ->default_val(42);
    cmd.add_option(
           "--checkpoint-step",
           "Save checkpoint every N iterations (omit to save only at end)")
        ->group("Runtime");
    cmd.add_option("--from-ckpt", "Run from checkpoint")
        ->group("Runtime")
        ->type_name("<CKPT>");
    cmd.add_option("-c,--chunk-size", "SNPs per chunk (controls memory usage)")
        ->group("Runtime")
        ->default_val(10000);
    cmd.add_option("-t,--threads", "CPU threads")
        ->group("Runtime")
        ->default_val(
            std::max(
                1, static_cast<int>(std::thread::hardware_concurrency() / 2)));
    cmd.add_flag(
           "--mmap", "Memory-mapped genotype I/O (lower RAM, may be slower)")
        ->group("Runtime");

    cmd.footer(
        "Docs:\n"
        "  https://gelex.readthedocs.io/en/latest/cli/mcmc.html");

    cmd.callback(
        [subcommand, &exit_code]()
        { exit_code = cli::execute_cli_command(*subcommand, mcmc_execute); });
}
