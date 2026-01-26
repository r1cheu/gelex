#include "fit_args.h"

#include <argparse.h>
#include <thread>

void setup_fit_args(argparse::ArgumentParser& cmd)
{
    cmd.add_description(
        "Fit genomic prediction models using Bayesian or GBLUP methods");

    // ================================================================
    // IO
    // ================================================================
    cmd.add_group("Data Files");
    cmd.add_argument("-p", "--pheno")
        .help("Phenotype file (TSV format: FID, IID, trait1, ...)")
        .metavar("<PHENOTYPE>")
        .required();
    cmd.add_argument("-b", "--bfile")
        .help("PLINK binary file prefix (.bed/.bim/.fam)")
        .metavar("<BFILE>")
        .required();
    cmd.add_argument("--qcovar")
        .default_value("")
        .help("Quantitative covariates (TSV: FID, IID, covar1, ...)");
    cmd.add_argument("--dcovar")
        .default_value("")
        .help("Discrete covariates (TSV: FID, IID, factor1, ...)");
    cmd.add_argument("-o", "--out")
        .help("Output file prefix")
        .metavar("<OUT>")
        .default_value("gelex");

    // ================================================================
    // Data Processing
    // ================================================================
    cmd.add_group("Processing Options");
    cmd.add_argument("--pheno-col")
        .help("Phenotype column index (0-based)")
        .default_value(2)
        .scan<'i', int>();
    cmd.add_argument("-c", "--chunk-size")
        .help("SNPs per chunk (controls memory usage)")
        .default_value(10000)
        .scan<'i', int>();
    cmd.add_argument("--iid-only")
        .help("Use only IID for sample matching (ignore FID)")
        .flag();

    // ================================================================
    // Model Configuration
    // ================================================================
    cmd.add_group("Model Configuration");
    cmd.add_argument("-m", "--method")
        .help(
            "Method: A/B/C/R/RR/GBLUP (+d for dominance, +pi to estimate "
            "mixture), e.g. RRd, Bdpi")
        .default_value("RR")
        .metavar("<METHOD>")
        .choices(
            "A",
            "Ad",
            "B",
            "Bpi",
            "Bd",
            "Bdpi",
            "C",
            "Cpi",
            "Cd",
            "Cdpi",
            "R",
            "Rd",
            "RR",
            "RRd")
        .required();

    cmd.add_argument("--scale")
        .help("Additive variance scales for BayesR (5 values)")
        .nargs(5)
        .default_value(std::vector<double>{0, 0.001, 0.01, 0.1, 1})
        .scan<'g', double>();
    cmd.add_argument("--pi")
        .help("Additive mixture proportions for BayesB/C/R")
        .nargs(2)
        .default_value(std::vector<double>{0.95, 0.05})
        .scan<'g', double>();
    cmd.add_argument("--dscale")
        .help("Dominance variance scales for BayesR (5 values)")
        .nargs(5)
        .default_value(std::vector<double>{0, 0.001, 0.01, 0.1, 1})
        .scan<'g', double>();
    cmd.add_argument("--dpi")
        .help("Dominance mixture proportions for BayesB/C/R")
        .nargs(2)
        .default_value(std::vector<double>{0.95, 0.05})
        .scan<'g', double>();

    // ================================================================
    // MCMC Parameters
    // ================================================================
    cmd.add_group("MCMC Configuration");
    cmd.add_argument("--iters")
        .help("Total MCMC iterations")
        .default_value(5000)
        .scan<'i', int>();
    cmd.add_argument("--burnin")
        .help("Burn-in iterations to discard")
        .default_value(4000)
        .scan<'i', int>();
    cmd.add_argument("--thin")
        .help("Thinning interval for samples")
        .default_value(1)
        .scan<'i', int>();
    cmd.add_argument("--chains")
        .help("Number of MCMC chains")
        .default_value(1)
        .scan<'i', int>();

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
    cmd.add_argument("--mmap")
        .help(
            "Use memory-mapped I/O for genotype matrix(much lower RAM, may be "
            "slower)")
        .flag();
}
