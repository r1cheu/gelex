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

#include "cli_helper.h"

void setup_simulate_args(argparse::ArgumentParser& cmd)
{
    cmd.add_description(
        "Simulate phenotypes based on genetic data and specified parameters");

    // ================================================================
    // Data Files
    // ================================================================
    cmd.add_group("Data Files");
    cmd.add_argument("-b", "--bfile")
        .help("PLINK binary file prefix (.bed/.bim/.fam)")
        .metavar("<BFILE>")
        .required();
    cmd.add_argument("-o", "--out")
        .help("Output file prefix for simulated phenotypes")
        .metavar("<OUT>")
        .default_value("sim.phen");

    // ================================================================
    // Simulation Parameters
    // ================================================================
    cmd.add_group("Simulation Parameters");
    cmd.add_argument("--h2")
        .help("Narrow-sense heritability (range: 0-1)")
        .default_value(0.5)
        .scan<'g', double>();
    cmd.add_argument("--d2")
        .help("Dominance variance proportion (range: 0-1, h2+d2<1)")
        .default_value(0.0)
        .scan<'g', double>();
    cmd.add_argument("--add-var")
        .help("Variances for additive effect classes")
        .metavar("<VARIANCES>")
        .nargs(argparse::nargs_pattern::at_least_one)
        .scan<'g', double>();
    cmd.add_argument("--add-prop")
        .help(
            "Proportions for additive effect classes "
            "(must match --add-var length, sum to 1)")
        .metavar("<PROPORTIONS>")
        .nargs(argparse::nargs_pattern::at_least_one)
        .scan<'g', double>();
    cmd.add_argument("--dom-var")
        .help("Variances for dominance effect classes")
        .metavar("<VARIANCES>")
        .nargs(argparse::nargs_pattern::at_least_one)
        .scan<'g', double>();
    cmd.add_argument("--dom-prop")
        .help(
            "Proportions for dominance effect classes "
            "(must match --dom-var length, sum to 1)")
        .metavar("<PROPORTIONS>")
        .nargs(argparse::nargs_pattern::at_least_one)
        .scan<'g', double>();
    cmd.add_argument("--seed")
        .help("Random seed for reproducibility (-1 for time-based)")
        .default_value(-1)
        .scan<'i', int>();

    cmd.add_epilog(
        gelex::cli::format_epilog(
            "{bg}Examples:{rs}\n"
            "  {gy}# Basic phenotype simulation{rs}\n"
            "  {bc}gelex simulate{rs} {cy}-b{rs} geno\n"
            "  {gy}# Custom heritability with dominance{rs}\n"
            "  {bc}gelex simulate{rs} {cy}-b{rs} geno "
            "{cy}--h2{rs} 0.3 {cy}--d2{rs} 0.1 {cy}--seed{rs} 42\n"
            "  {gy}# Mixture normal effect sizes (BayesR-style){rs}\n"
            "  {bc}gelex simulate{rs} {cy}-b{rs} geno "
            "{cy}--add-var{rs} 0 0.01 0.001 0.0001 "
            "{cy}--add-prop{rs} 0.90 0.003 0.103 0.894 "
            "{cy}--h2{rs} 0.5 {cy}--seed{rs} 42"));
}
