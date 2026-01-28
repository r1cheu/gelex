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

#include "simulation_args.h"

#include <argparse.h>

void setup_simulation_args(argparse::ArgumentParser& cmd)
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
    cmd.add_argument("--causal")
        .help(
            "Causal variants file (TSV: SNP ID per line, optional effect size)")
        .metavar("<CAUSAL>")
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
    cmd.add_argument("--seed")
        .help("Random seed for reproducibility (-1 for time-based)")
        .default_value(-1)
        .scan<'i', int>();
}
