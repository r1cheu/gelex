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

#include "assoc_args.h"

#include <argparse.h>
#include <thread>

#include "cli/cli_helper.h"

auto setup_assoc_args(argparse::ArgumentParser& cmd) -> void
{
    cmd.add_description(
        "Perform genome-wide association study using mixed linear model");

    // ================================================================
    // IO
    // ================================================================
    cmd.add_group("Data Files");
    cmd.add_argument("-p", "--pheno")
        .help("Phenotype file (TSV format: FID, IID, trait1, ...)")
        .metavar("<PHENOTYPE>")
        .required();
    cmd.add_argument("--pheno-col")
        .help("Phenotype column index (0-based)")
        .default_value(2)
        .scan<'i', int>();
    cmd.add_argument("-b", "--bfile")
        .help("PLINK binary file prefix (.bed/.bim/.fam)")
        .metavar("<BFILE>")
        .required();
    cmd.add_argument("--grm")
        .help("GRM file prefix(es). Can specify multiple GRMs.")
        .metavar("<GRM>")
        .nargs(argparse::nargs_pattern::at_least_one)
        .required();
    cmd.add_argument("--qcovar")
        .help("Quantitative covariates (TSV: FID, IID, covar1, ...)");
    cmd.add_argument("--dcovar")
        .help("Discrete covariates (TSV: FID, IID, factor1, ...)");
    cmd.add_argument("-o", "--out")
        .help("Output file prefix")
        .metavar("<OUT>")
        .default_value("gelex");
    // ================================================================
    // REML Configuration
    // ================================================================
    cmd.add_group("REML Options");
    cmd.add_argument("--max-iter")
        .help("Max iteration in REML process")
        .default_value(100)
        .scan<'i', int>();
    cmd.add_argument("--tol")
        .help("tolerance for convergence in REML process")
        .default_value(1e-6)
        .scan<'g', double>();

    // ================================================================
    // Data Processing
    // ================================================================
    cmd.add_group("Processing Options");
    cmd.add_argument("-c", "--chunk-size")
        .help("SNPs per chunk for association testing")
        .default_value(10000)
        .scan<'i', int>();
    cmd.add_argument("--loco")
        .help("Enable Leave-One-Chromosome-Out (LOCO) mode")
        .flag();

    // ================================================================
    // Model Configuration
    // ================================================================
    cmd.add_group("Model Configuration");
    cmd.add_argument("--model")
        .help(
            "Association model: a for additive association test, d for "
            "dominance association test")
        .default_value("a")
        .metavar("<MODEL>")
        .choices("a", "d");
    cmd.add_argument("--geno-method", "--gm")
        .help(
            "Genotype method: StandardizeHWE(SH), CenterHWE(CH),"
            " OrthStandardizeHWE(OSH), OrthCenterHWE(OCH),"
            " Standardize(S), Center(C), OrthStandardize(OS), OrthCenter(OC)")
        .default_value(std::string("OCH"))
        .metavar("<STR>");
    cmd.add_argument("--transform")
        .help(
            "Phenotype transformation: none (default), dint (Direct INT), "
            "iint (Indirect INT)")
        .default_value("none")
        .metavar("<TRANSFORM>")
        .choices("none", "dint", "iint");
    cmd.add_argument("--int-offset")
        .help("INT offset parameter k (default: 3/8 Blom offset)")
        .default_value(3.0 / 8.0)
        .scan<'g', double>();
    // ================================================================
    // Performance
    // ================================================================
    cmd.add_group("Performance");
    cmd.add_argument("-t", "--threads")
        .help("Number of CPU threads to use")
        .default_value(
            std::max(
                1, static_cast<int>(std::thread::hardware_concurrency() / 2)))
        .scan<'i', int>();

    cmd.add_epilog(
        gelex::cli::format_epilog(
            "{bg}Example:{rs}\n"
            "  {bc}gelex assoc{rs} {cy}-p{rs} pheno.tsv {cy}-b{rs} geno "
            "{cy}--grm{rs} grm_prefix\n\n"
            "{bg}Docs:{rs}\n"
            "  https://gelex.readthedocs.io/en/latest/cli/assoc.html"));
}
