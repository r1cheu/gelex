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

#include "predict_command.h"

#include <argparse.h>

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
