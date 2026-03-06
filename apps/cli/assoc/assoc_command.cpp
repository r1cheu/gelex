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

#include "assoc_command.h"

#include <argparse.h>

#include "assoc_config.h"
#include "assoc_reporter.h"
#include "cli/cli_helper.h"
#include "cli/data_pipe_config.h"
#include "cli/data_pipe_reporter.h"
#include "cli/reml_reporter.h"
#include "gelex/infra/logging/assoc_event.h"
#include "gelex/pipeline/assoc_loco_engine.h"
#include "gelex/pipeline/assoc_normal_engine.h"
#include "gelex/pipeline/grm_pipe.h"
#include "gelex/pipeline/pheno_pipe.h"
#include "gelex/types/genetic_effect_type.h"

auto assoc_execute(argparse::ArgumentParser& cmd) -> int
{
    bool loco = cmd.get<bool>("--loco");
    int threads = cmd.get<int>("--threads");
    gelex::cli::setup_parallelization(threads);

    auto pheno_config = gelex::cli::make_pheno_config(cmd);
    pheno_config.transform_type
        = gelex::cli::parse_transform_type(cmd.get("--transform"));
    pheno_config.int_offset = cmd.get<double>("--int-offset");

    gelex::cli::AssocReporter reporter;
    gelex::cli::DataPipeReporter data_reporter;

    reporter.on_event(
        gelex::AssocConfigLoadedEvent{
            .model_type = cmd.get("--model") == "a" ? gelex::ModelType::A
                                                    : gelex::ModelType::D,
            .loco = loco,

            .geno_method = gelex::cli::parse_genotype_process_method(
                cmd.get<std::string>("--geno-method")),

            .max_iter = cmd.get<int>("--max-iter"),
            .tol = cmd.get<double>("--tol"),
        });

    std::vector<std::filesystem::path> grm_paths;
    for (const auto& p : cmd.get<std::vector<std::string>>("--grm"))
    {
        grm_paths.emplace_back(p);
    }

    gelex::PhenoPipe pheno(pheno_config, data_reporter.as_observer());
    gelex::GrmPipe grm(grm_paths, data_reporter.as_observer());
    pheno.load(grm.sample_id_sets());
    grm.load(pheno.sample_manager());

    auto config = gelex::cli::make_assoc_config(cmd);

    if (loco)
    {
        gelex::AssocLocoEngine engine(std::move(config));
        engine.run(pheno, grm, reporter.as_observer());
    }
    else
    {
        gelex::cli::RemlReporter reml_reporter;
        gelex::AssocNormalEngine engine(std::move(config));
        engine.run(
            pheno, grm, reporter.as_observer(), reml_reporter.as_observer());
    }

    return 0;
}
