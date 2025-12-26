#include "gelex/cli/predict_command.h"

#include <thread>

#include "gelex/data/bed_pipe.h"
#include "gelex/logger.h"
#include "predict/predict_engine.h"

void predict_command(argparse::ArgumentParser& cmd)
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
    cmd.add_argument("-c", "--covar-eff")
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
    cmd.add_argument("--chunk-size")
        .help("SNPs per chunk (controls memory usage)")
        .default_value(10000)
        .scan<'i', int>();
}

int predict_execute(argparse::ArgumentParser& predict)
{
    auto logger = gelex::logging::get();

    gelex::PredictEngine::Config config;
    config.bed_path = gelex::BedPipe::format_bed_path(predict.get("bfile"));
    config.snp_effect_path = predict.get("--snp-eff");
    config.covar_effect_path = predict.get("--covar-eff");
    config.qcovar_path = predict.get("--qcovar");
    config.dcovar_path = predict.get("--dcovar");
    config.output_path = predict.get("--out");
    config.iid_only = predict.get<bool>("--iid-only");

    try
    {
        config.validate();
    }
    catch (const std::exception& e)
    {
        if (logger)
        {
            logger->error("Configuration validation failed: {}", e.what());
        }
        else
        {
            std::cerr << "[error] " << e.what() << "\n";
        }
        return 1;
    }

    try
    {
        gelex::PredictEngine engine(config);
        engine.run();
    }
    catch (const std::exception& e)
    {
        if (logger)
        {
            logger->error("Prediction failed: {}", e.what());
        }
        else
        {
            std::cerr << "[error] " << e.what() << "\n";
        }
        return 1;
    }

    if (logger)
    {
        logger->info("Prediction completed successfully");
    }

    return 0;
}
