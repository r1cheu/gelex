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

#include "gelex/pipeline/predict_engine.h"
#include "predict_config.h"
#include "predict_reporter.h"

auto predict_execute(argparse::ArgumentParser& predict) -> int
{
    auto configs = gelex::cli::make_predict_configs(predict);
    gelex::cli::PredictReporter reporter;

    gelex::PredictParamsPipe params_pipe(
        configs.params, reporter.as_observer());

    gelex::PredictDataPipe data_pipe(
        configs.data,
        params_pipe.snp_effects(),
        params_pipe.covar_params(),
        reporter.as_observer());

    gelex::PredictEngine engine(std::move(configs.engine));
    engine.run(
        std::move(params_pipe), std::move(data_pipe), reporter.as_observer());

    return 0;
}
