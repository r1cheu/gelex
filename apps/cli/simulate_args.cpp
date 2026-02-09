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
        .help("Variances for additive effect classes (default: 0.01)")
        .metavar("<VARIANCES>")
        .nargs(argparse::nargs_pattern::at_least_one)
        .default_value(std::vector<double>{0.01})
        .scan<'g', double>();
    cmd.add_argument("--add-prop")
        .help(
            "Proportions for additive effect classes "
            "(must match --add-var length, sum to 1, default: 1.0)")
        .metavar("<PROPORTIONS>")
        .nargs(argparse::nargs_pattern::at_least_one)
        .default_value(std::vector<double>{1.0})
        .scan<'g', double>();
    cmd.add_argument("--dom-var")
        .help("Variances for dominance effect classes (default: 0.01)")
        .metavar("<VARIANCES>")
        .nargs(argparse::nargs_pattern::at_least_one)
        .default_value(std::vector<double>{0.01})
        .scan<'g', double>();
    cmd.add_argument("--dom-prop")
        .help(
            "Proportions for dominance effect classes "
            "(must match --dom-var length, sum to 1, default: 1.0)")
        .metavar("<PROPORTIONS>")
        .nargs(argparse::nargs_pattern::at_least_one)
        .default_value(std::vector<double>{1.0})
        .scan<'g', double>();
    cmd.add_argument("--intercept")
        .help("Intercept (mean) term added to phenotypes")
        .default_value(0.0)
        .scan<'g', double>();
    cmd.add_argument("--seed")
        .help("Random seed for reproducibility")
        .default_value(42)
        .scan<'i', int>();

    cmd.add_epilog(
        gelex::cli::format_epilog(
            "{bg}Example:{rs}\n"
            "  {bc}gelex simulate{rs} {cy}-b{rs} geno\n\n"
            "{bg}Docs:{rs}\n"
            "  https://gelex.readthedocs.io/en/latest/cli/simulate.html"));
}
