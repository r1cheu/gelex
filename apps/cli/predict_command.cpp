#include "predict_command.h"

#include <thread>

#include "gelex/data/bed_pipe.h"
#include "gelex/logger.h"
#include "predict/predict_engine.h"

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
