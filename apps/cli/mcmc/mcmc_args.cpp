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

#include "mcmc_args.h"

#include <algorithm>
#include <string>
#include <thread>

#include <argparse.h>

#include "cli/cli_helper.h"

auto setup_mcmc_args(argparse::ArgumentParser& cmd) -> void
{
    cmd.add_description(
        "Fit genomic prediction models using MCMC (Gibbs sampling)");

    cmd.add_group("I/O");
    cmd.add_argument("-p", "--pheno")
        .help("Phenotype file (TSV format: FID, IID, trait1, ...)")
        .metavar("<PHENOTYPE>")
        .required();
    auto& geno_group = cmd.add_mutually_exclusive_group(true);
    geno_group.add_argument("-b", "--bfile")
        .help("PLINK binary file prefix (.bed/.bim/.fam)")
        .metavar("<BFILE>");
    geno_group.add_argument("-g", "--gfile")
        .help("Gelex binary genotype file prefix")
        .metavar("<GFILE>");
    cmd.add_argument("--qcovar")
        .help("Quantitative covariates (TSV: FID, IID, covar1, ...)");
    cmd.add_argument("--dcovar")
        .help("Discrete covariates (TSV: FID, IID, factor1, ...)");
    cmd.add_argument("-o", "--out")
        .help("Output file prefix")
        .metavar("<OUT>")
        .default_value("gelex");

    cmd.add_group("Processing");
    cmd.add_argument("--pheno-col")
        .help("Phenotype column index (0-based)")
        .default_value(2)
        .scan<'i', int>();
    cmd.add_argument("-c", "--chunk-size")
        .help("SNPs per chunk (controls memory usage)")
        .default_value(10000)
        .scan<'i', int>();
    cmd.add_argument("--geno-method", "--gm")
        .help(
            "Genotype method: StandardizeHWE(SH), CenterHWE(CH),"
            " OrthStandardizeHWE(OSH), OrthCenterHWE(OCH),"
            " Standardize(S), Center(C), OrthStandardize(OS), OrthCenter(OC),"
            " NOIAStandardize(NS), NOIACenter(NC)")
        .default_value(std::string("OSH"))
        .metavar("<STR>");

    cmd.add_group("Model");
    cmd.add_argument("-m", "--method")
        .help("Bayesian method: RR, A, B, C, R, CD")
        .default_value("RR")
        .metavar("<METHOD>")
        .choices("A", "B", "C", "R", "RR", "CD")
        .required();
    cmd.add_argument("--mode")
        .help(
            "Genetic effect mode: A (additive only), D (dominance only), "
            "AD (additive + dominance)")
        .default_value(std::string("A"))
        .choices("A", "D", "AD")
        .metavar("<MODE>");
    cmd.add_argument("--random-variance-proportion")
        .help("Variance fraction for non-SNP random effects (0, 1)")
        .scan<'g', double>();
    cmd.add_argument("--additive-h2")
        .help("Additive heritability (0, 1)")
        .scan<'g', double>();
    cmd.add_argument("--dominance-h2")
        .help("Dominance heritability (0, 1)")
        .scan<'g', double>();
    cmd.add_argument("--additive-pi")
        .help("Additive mixture proportions (B/C/R)")
        .nargs(argparse::nargs_pattern::at_least_one)
        .scan<'g', double>();
    cmd.add_argument("--dominance-pi")
        .help("Dominance mixture proportions (B/C/R)")
        .nargs(argparse::nargs_pattern::at_least_one)
        .scan<'g', double>();
    cmd.add_argument("--additive-multiplier")
        .help("Additive variance multipliers (R)")
        .nargs(argparse::nargs_pattern::at_least_one)
        .scan<'g', double>();
    cmd.add_argument("--dominance-multiplier")
        .help("Dominance variance multipliers (R)")
        .nargs(argparse::nargs_pattern::at_least_one)
        .scan<'g', double>();
    cmd.add_argument("--joint-pi")
        .help("Joint allocation proportions (CD)")
        .nargs(argparse::nargs_pattern::at_least_one)
        .scan<'g', double>();
    cmd.add_argument("--estimate-additive-pi")
        .help("Sample additive mixture proportions")
        .flag();
    cmd.add_argument("--estimate-dominance-pi")
        .help("Sample dominance mixture proportions")
        .flag();
    cmd.add_argument("--estimate-joint-pi")
        .help("Sample joint allocation proportions (CD)")
        .flag();

    cmd.add_group("MCMC");
    cmd.add_argument("--iters")
        .help("MCMC iterations to run")
        .default_value(5000)
        .scan<'i', int>();
    cmd.add_argument("--burn-in")
        .help("Burn-in iterations to discard")
        .default_value(3000)
        .scan<'i', int>();
    cmd.add_argument("--thin")
        .help("Thinning interval for samples")
        .default_value(1)
        .scan<'i', int>();
    cmd.add_argument("--seed")
        .help("Random seed for MCMC")
        .default_value(42)
        .scan<'i', int>();
    cmd.add_argument("--checkpoint-step")
        .help("Save checkpoint every N iterations (omit to save only at end)")
        .scan<'i', int>();
    cmd.add_argument("--from-ckpt")
        .help("Run from checkpoint")
        .metavar("<CKPT>");

    cmd.add_group("Performance");
    cmd.add_argument("-t", "--threads")
        .help("CPU threads")
        .default_value(
            std::max(
                1, static_cast<int>(std::thread::hardware_concurrency() / 2)))
        .scan<'i', int>();
    cmd.add_argument("--mmap")
        .help("Memory-mapped genotype I/O (lower RAM, may be slower)")
        .flag();

    cmd.add_epilog(
        gelex::cli::format_epilog(
            "{bg}Example:{rs}\n"
            "  {bc}gelex mcmc{rs} {cy}-p{rs} pheno.tsv {cy}-b{rs} geno "
            "{cy}-m{rs} R {cy}--mode AD{rs}\n\n"
            "{bg}Docs:{rs}\n"
            "  https://gelex.readthedocs.io/en/latest/cli/mcmc.html"));
}
