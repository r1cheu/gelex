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
#include <filesystem>
#include <string>
#include <vector>

#include "cli/cli_helper.h"
#include "cli/data_pipe_config.h"
#include "cli/data_pipe_reporter.h"
#include "fit_config.h"
#include "fit_reporter.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/genotype/bed_path.h"
#include "gelex/data/reader.h"
#include "gelex/infra/logging/data_pipe_event.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/pipeline/fit_engine.h"
#include "gelex/pipeline/geno_pipe.h"
#include "gelex/pipeline/pheno_pipe.h"
#include "gelex/types/genotype_process_method.h"

auto fit_execute(argparse::ArgumentParser& fit) -> int
{
    auto fit_config = gelex::cli::make_fit_config(fit);
    auto [pheno_config, geno_config]
        = gelex::cli::make_fit_data_configs(fit, fit.get<bool>("--mmap"));

    auto model_type = fit_config.method.dominance ? gelex::GeneticMode::AD
                                                  : gelex::GeneticMode::A;

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
            .n_burn_in = static_cast<int>(fit_config.mcmc_params.n_burn_in),
            .seed = fit_config.seed,
        });

    auto bed_path = gelex::format_bed_path(fit.get<std::string>("--bfile"));
    auto fam_index
        = gelex::read_fam(
              std::filesystem::path(bed_path).replace_extension(".fam"))
              .index();

    gelex::PhenoPipe pheno(pheno_config, data_reporter.as_observer());
    pheno.load();

    std::vector<const gelex::df::Index<std::string>*> all_indices{
        &fam_index, &pheno.pheno_index()};
    all_indices.append_range(pheno.covar_indices());
    auto common = gelex::df::intersect<std::string>(all_indices);

    gelex::notify(
        data_reporter.as_observer(),
        gelex::IntersectionEvent{.common_samples = common.size()});

    pheno.gather(common);

    gelex::GenoPipe geno(geno_config, data_reporter.as_observer());
    geno.load(common);

    gelex::FitEngine engine(std::move(fit_config));
    engine.run(std::move(pheno), std::move(geno), reporter.as_observer());

    return 0;
}
