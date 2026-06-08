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

#include <argparse.h>

#include "cli/cli_helper.h"

auto setup_reml_args(argparse::ArgumentParser& cmd) -> void
{
    cmd.add_description(
        "Perform REML variance component estimation using average information "
        "algorithm");

    cmd.add_group("I/O");
    cmd.add_argument("-p", "--pheno")
        .help("Phenotype file (TSV format: FID, IID, trait1, ...)")
        .metavar("<PHENOTYPE>")
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
    cmd.add_argument("--rand").help(
        "Random effects (TSV: FID, IID, effect1, ...)");
    cmd.add_argument("-o", "--out")
        .help("Output file prefix")
        .metavar("<OUT>")
        .default_value("gelex");

    cmd.add_group("Phenotype");
    cmd.add_argument("--pheno-col")
        .help("Phenotype column index (0-based)")
        .default_value(2)
        .scan<'i', int>();
    cmd.add_argument("--geno-method", "--gm")
        .help(
            "Genotype processing: SH, CH, OSH, OCH, S, C, OS, OC, NS, NC "
            "(prefix: O=orth, N=NOIA; suffix: H=HWE-based)")
        .default_value(std::string("OCH"))
        .metavar("<STR>");

    cmd.add_group("REML");
    cmd.add_argument("--max-iter")
        .help("Max iterations")
        .default_value(100)
        .scan<'i', int>();
    cmd.add_argument("--tol")
        .help("Convergence tolerance")
        .default_value(1e-6)
        .scan<'g', double>();

    cmd.add_group("Performance");
    cmd.add_argument("-t", "--threads")
        .help("CPU threads")
        .default_value(
            std::max(
                1, static_cast<int>(std::thread::hardware_concurrency() / 2)))
        .scan<'i', int>();

    cmd.add_epilog(
        gelex::cli::format_epilog(
            "{bg}Example:{rs}\n"
            "  {bc}gelex reml{rs} {cy}-p{rs} pheno.tsv geno {cy}--grm{rs} "
            "grm_prefix\n\n"
            "{bg}Docs:{rs}\n"
            "  https://gelex.readthedocs.io/en/latest/cli/reml.html"));
}
