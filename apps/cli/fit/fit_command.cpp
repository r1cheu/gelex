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
#include "cli/config_factory.h"
#include "fit_config.h"
#include "fit_reporter.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/pipeline/fit_engine.h"

int fit_execute(argparse::ArgumentParser& fit)
{
    auto config = gelex::cli::make_config<FitConfig>(fit);
    gelex::cli::FitReporter reporter;

    gelex::cli::setup_parallelization(config.threads);

    reporter.on_event(
        gelex::FitConfigLoadedEvent{
            .method = config.method,
            .use_dominance = config.use_dominance,
            .n_iters = static_cast<int>(config.mcmc_params.n_iters),
            .n_burnin = static_cast<int>(config.mcmc_params.n_burnin),
            .seed = config.seed,
            .threads = config.threads,
        });

    gelex::FitEngine engine({
        .bed_path = config.bed_path,
        .out_prefix = config.out_prefix,
        .method = config.method,
        .use_dominance = config.use_dominance,
        .seed = config.seed,
        .mcmc_params = config.mcmc_params,
        .phenotype_path = config.phenotype_path,
        .phenotype_column = config.phenotype_column,
        .use_mmap = config.use_mmap,
        .chunk_size = config.chunk_size,
        .qcovar_path = config.qcovar_path,
        .dcovar_path = config.dcovar_path,
        .genotype_method = config.genotype_method,
        .pi = config.pi,
        .dpi = config.dpi,
        .scale = config.scale,
        .dscale = config.dscale,
    });

    engine.run(
        [&reporter](const gelex::FitEvent& event)
        {
            std::visit(
                [&reporter](const auto& e) { reporter.on_event(e); }, event);
        });

    return 0;
}
