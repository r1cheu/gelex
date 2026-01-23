#include "grm_args.h"

#include <thread>

#include <argparse.h>

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
    cmd.add_argument("--method")
        .help("GRM computation method: su, yang, zeng, vitezica")
        .metavar("<METHOD>")
        .default_value(std::string("su"));
    cmd.add_argument("--chunk")
        .help("Chunk size for memory-efficient computation")
        .metavar("<SIZE>")
        .default_value(10000)
        .scan<'i', int>();
    cmd.add_argument("--threads")
        .help("Number of threads (-1 for all cores)")
        .metavar("<N>")
        .default_value(std::thread::hardware_concurrency() / 2)
        .scan<'i', int>();
    cmd.add_argument("--additive").help("Compute additive GRM").flag();
    cmd.add_argument("--dominant").help("Compute dominance GRM").flag();
    cmd.add_argument("--loco").help("Compute GRM for each chromosome").flag();
}
