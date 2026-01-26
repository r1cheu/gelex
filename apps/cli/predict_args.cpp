#include "predict_args.h"

#include <argparse.h>

void setup_predict_args(argparse::ArgumentParser& cmd)
{
    cmd.add_description(
        "Generate genomic predictions using fitted SNP effects");

    // ================================================================
    // Data Files
    // ================================================================
    cmd.add_group("Data Files");
    cmd.add_argument("-b", "--bfile")
        .help("PLINK binary file prefix for prediction data (.bed/.bim/.fam)")
        .metavar("<BFILE>")
        .required();
    cmd.add_argument("-e", "--snp-eff")
        .help("SNP effects file (.snp.eff)")
        .metavar("<SNP_EFF>")
        .required();
    cmd.add_argument("--covar-eff")
        .help("Covariate effects file (.param)")
        .metavar("<COVAR_EFF>");
    cmd.add_argument("--qcovar")
        .help("Quantitative covariates file (TSV: FID, IID, covar1, ...)")
        .default_value("")
        .metavar("<QCOVAR>");
    cmd.add_argument("--dcovar")
        .help("Discrete covariates file (TSV: FID, IID, factor1, ...)")
        .default_value("")
        .metavar("<DCOVAR>");
    cmd.add_argument("-o", "--out")
        .help("Output file path for predictions")
        .metavar("<OUT>")
        .required();

    // ================================================================
    // Processing Options
    // ================================================================
    cmd.add_group("Processing Options");
    cmd.add_argument("--iid-only")
        .help("Use only IID for sample matching (ignore FID)")
        .flag();
    cmd.add_argument("-c", "--chunk-size")
        .help("SNPs per chunk (controls memory usage)")
        .default_value(10000)
        .scan<'i', int>();
}
