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

#include "fit_command.h"

#include <argparse.h>
#include <variant>

#include "cli/cli_helper.h"
#include "cli/data_pipe_config.h"
#include "cli/data_pipe_reporter.h"
#include "fit_config.h"
#include "fit_reporter.h"
#include "gelex/infra/logging/data_pipe_event.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/pipeline/data_pipe.h"
#include "gelex/pipeline/fit_engine.h"

auto fit_execute(argparse::ArgumentParser& fit) -> int
{
    auto fit_config = gelex::cli::make_fit_config(fit);
    auto data_config
        = gelex::cli::make_data_config(fit, fit.get<bool>("--mmap"));

    data_config.model_type = gelex::cli::has_dominance(fit_config.method)
                                 ? gelex::ModelType::AD
                                 : gelex::ModelType::A;

    int threads = fit.get<int>("--threads");
    gelex::cli::FitReporter reporter;
    gelex::cli::DataPipeReporter data_reporter;
    gelex::cli::setup_parallelization(threads);

    reporter.on_event(
        gelex::FitConfigLoadedEvent{
            .method = fit_config.method,
            .n_iters = static_cast<int>(fit_config.mcmc_params.n_iters),
            .n_burnin = static_cast<int>(fit_config.mcmc_params.n_burnin),
            .seed = fit_config.seed,
            .threads = threads,
        });

    gelex::DataPipe data(
        data_config,
        [&data_reporter](const gelex::DataPipeEvent& e)
        {
            std::visit([&](const auto& ev) { data_reporter.on_event(ev); }, e);
        });

    data.load();

    gelex::FitEngine engine(std::move(fit_config));

    engine.run(
        std::move(data),
        [&reporter](const gelex::FitEvent& e)
        { std::visit([&](const auto& ev) { reporter.on_event(ev); }, e); });

    return 0;
}
