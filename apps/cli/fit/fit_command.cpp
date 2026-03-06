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

#include "cli/cli_helper.h"
#include "cli/data_pipe_config.h"
#include "cli/data_pipe_reporter.h"
#include "fit_config.h"
#include "fit_reporter.h"
#include "gelex/infra/logging/data_pipe_event.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/pipeline/fit_engine.h"
#include "gelex/pipeline/geno_pipe.h"
#include "gelex/pipeline/pheno_pipe.h"

auto fit_execute(argparse::ArgumentParser& fit) -> int
{
    auto fit_config = gelex::cli::make_fit_config(fit);
    auto [pheno_config, geno_config]
        = gelex::cli::make_fit_data_configs(fit, fit.get<bool>("--mmap"));

    auto model_type = gelex::cli::has_dominance(fit_config.method)
                          ? gelex::ModelType::AD
                          : gelex::ModelType::A;

    geno_config.model_type = model_type;

    int threads = fit.get<int>("--threads");
    gelex::cli::FitReporter reporter;
    gelex::cli::DataPipeReporter data_reporter;
    gelex::cli::setup_parallelization(threads);

    reporter.on_event(
        gelex::FitConfigLoadedEvent{
            .method = fit_config.method,
            .model_type = model_type,
            .n_iters = static_cast<int>(fit_config.mcmc_params.n_iters),
            .n_burnin = static_cast<int>(fit_config.mcmc_params.n_burnin),
            .seed = fit_config.seed,
        });

    gelex::PhenoPipe pheno(pheno_config, data_reporter.as_observer());
    pheno.load();

    gelex::GenoPipe geno(geno_config, data_reporter.as_observer());
    geno.load(pheno.sample_manager());

    gelex::FitEngine engine(std::move(fit_config));

    engine.run(std::move(pheno), std::move(geno), reporter.as_observer());

    return 0;
}
