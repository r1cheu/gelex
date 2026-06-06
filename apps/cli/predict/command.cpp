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

#include "command.h"

#include <argparse.h>
#include <utility>

#include "config.h"
#include "gelex/engine/predict.h"
#include "gelex/infra/logging/predict_event.h"
#include "reporter.h"

auto predict_execute(argparse::ArgumentParser& predict) -> int
{
    auto config = gelex::cli::make_predict_config(predict);
    gelex::cli::PredictReporter reporter;
    reporter.on_event(gelex::PredictBannerEvent{});
    gelex::PredictEngine engine(std::move(config));
    engine.run(reporter.as_observer());
    return 0;
}
