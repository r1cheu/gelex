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

#include "predict_args.h"

#include <argparse.h>

#include "cli_helper.h"

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

    cmd.add_epilog(
        gelex::cli::format_epilog(
            "{bg}Examples:{rs}\n"
            "  {gy}# Basic genomic prediction{rs}\n"
            "  {bc}gelex predict{rs} {cy}-b{rs} geno {cy}-e{rs} model.snp.eff "
            "{cy}-o{rs} pred.tsv\n"
            "  {gy}# Prediction with covariate effects{rs}\n"
            "  {bc}gelex predict{rs} {cy}-b{rs} geno {cy}-e{rs} model.snp.eff "
            "{cy}--covar-eff{rs} model.param {cy}--qcovar{rs} age.txt "
            "{cy}-o{rs} pred.tsv"));
}
