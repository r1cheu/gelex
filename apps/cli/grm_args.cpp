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

#include "grm_args.h"

#include <thread>

#include <argparse.h>

#include "cli_helper.h"

void setup_grm_args(argparse::ArgumentParser& cmd)
{
    cmd.add_description(
        "Compute genomic relationship matrix (GRM) from PLINK "
        "binary files and output in GCTA format");

    // ================================================================
    // Data Files
    // ================================================================
    cmd.add_group("Data Files");
    cmd.add_argument("-b", "--bfile")
        .help("PLINK binary file prefix (.bed/.bim/.fam)")
        .metavar("<BFILE>")
        .required();
    cmd.add_argument("-o", "--out")
        .help("Output file prefix")
        .metavar("<OUT>")
        .default_value(std::string("grm"));

    // ================================================================
    // GRM Options
    // ================================================================
    cmd.add_group("GRM Options");
    cmd.add_argument("-m", "--method")
        .help("GRM computation method: su, yang, zeng, vitezica")
        .metavar("<METHOD>")
        .default_value(std::string("yang"));
    cmd.add_argument("-c", "--chunk-size")
        .help("Chunk size for memory-efficient computation")
        .metavar("<SIZE>")
        .default_value(10000)
        .scan<'i', int>();
    cmd.add_argument("-t", "--threads")
        .help("Number of threads (-1 for all cores)")
        .metavar("<N>")
        .default_value(
            static_cast<int>(std::thread::hardware_concurrency() / 2))
        .scan<'i', int>();
    cmd.add_argument("--add").help("Compute additive GRM").flag();
    cmd.add_argument("--dom").help("Compute dominance GRM").flag();
    cmd.add_argument("--loco").help("Compute GRM for each chromosome").flag();

    cmd.add_epilog(
        gelex::cli::format_epilog(
            "{bg}Examples:{rs}\n"
            "  {gy}# Compute additive GRM{rs}\n"
            "  {bc}gelex grm{rs} {cy}-b{rs} geno {cy}--add{rs}\n\n"
            "  {gy}# Compute dominance GRM with custom output{rs}\n"
            "  {bc}gelex grm{rs} {cy}-b{rs} geno {cy}--dom{rs} {cy}-o{rs} "
            "dom_grm\n\n"
            "  {gy}# LOCO GRM (one per chromosome){rs}\n"
            "  {bc}gelex grm{rs} {cy}-b{rs} geno {cy}--add{rs} "
            "{cy}--loco{rs}\n\n"
            "  {gy}# Compute both additive and dominance GRMs{rs}\n"
            "  {bc}gelex grm{rs} {cy}-b{rs} geno {cy}--add{rs} "
            "{cy}--dom{rs}\n\n"
            "  {gy}# Use specific method and threads{rs}\n"
            "  {bc}gelex grm{rs} {cy}-b{rs} geno {cy}--add{rs} {cy}-m{rs} "
            "vitezica {cy}-t{rs} 8"));
}
