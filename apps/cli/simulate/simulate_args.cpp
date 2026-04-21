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

#include "simulate_args.h"

#include <argparse.h>

#include "cli/cli_helper.h"

auto setup_simulate_args(argparse::ArgumentParser& cmd) -> void
{
    cmd.add_description(
        "Simulate phenotypes based on genetic data and specified parameters");

    cmd.add_group("I/O");
    cmd.add_argument("-b", "--bfile")
        .help("PLINK binary file prefix (.bed/.bim/.fam)")
        .metavar("<BFILE>")
        .required();
    cmd.add_argument("-o", "--out")
        .help("Output file prefix")
        .metavar("<OUT>")
        .default_value("sim.phen");

    cmd.add_group("Additive");
    cmd.add_argument("--h2")
        .help("Narrow-sense heritability (0, 1); omit to disable additive")
        .scan<'g', double>();
    cmd.add_argument("--add-var")
        .help("Variances for additive effect classes")
        .metavar("<VARIANCES>")
        .nargs(argparse::nargs_pattern::at_least_one)
        .scan<'g', double>();
    cmd.add_argument("--add-n")
        .help(
            "SNP counts for additive effect classes (must match --add-var "
            "length)")
        .metavar("<COUNTS>")
        .nargs(argparse::nargs_pattern::at_least_one)
        .scan<'i', int>();

    cmd.add_group("Dominance");
    cmd.add_argument("--d2")
        .help("Dominance variance proportion (0, 1); h2+d2<1 in AD mode")
        .scan<'g', double>();
    cmd.add_argument("--dom-var")
        .help("Variances for dominance effect classes")
        .metavar("<VARIANCES>")
        .nargs(argparse::nargs_pattern::at_least_one)
        .scan<'g', double>();
    cmd.add_argument("--dom-n")
        .help(
            "SNP counts for dominance effect classes (must match --dom-var "
            "length)")
        .metavar("<COUNTS>")
        .nargs(argparse::nargs_pattern::at_least_one)
        .scan<'i', int>();
    cmd.add_argument("--dom-pos-prob")
        .help(
            "Probability of positive dominance effects; enables "
            "truncated-normal sampling [0, 1]")
        .metavar("<PROB>")
        .scan<'g', double>();

    cmd.add_group("Model");
    cmd.add_argument("--geno-method", "--gm")
        .help(
            "Genotype processing: SH, CH, OSH, OCH, S, C, OS, OC, NS, NC "
            "(prefix: O=orth, N=NOIA; suffix: H=HWE-based)")
        .default_value(std::string("OS"))
        .metavar("<STR>");
    cmd.add_argument("--seed")
        .help("Random seed")
        .default_value(42)
        .scan<'i', int>();

    cmd.add_epilog(
        gelex::cli::format_epilog(
            "{bg}Example:{rs}\n"
            "  {bc}gelex simulate{rs} {cy}-b{rs} geno\n\n"
            "{bg}Docs:{rs}\n"
            "  https://gelex.readthedocs.io/en/latest/cli/simulate.html"));
}
