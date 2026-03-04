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
#include <variant>

#include "assoc_config.h"
#include "cli/cli_helper.h"
#include "cli/data_pipe_config.h"
#include "cli/data_pipe_reporter.h"
#include "gelex/data/genotype/bed_pipe.h"
#include "gelex/data/loader/bim_loader.h"
#include "gelex/infra/logging/data_pipe_event.h"
#include "gelex/pipeline/data_pipe.h"
#include "gwas_runner.h"

auto assoc_execute(argparse::ArgumentParser& cmd) -> int
{
    auto assoc_config = gelex::cli::make_assoc_config(cmd);
    auto data_config = gelex::cli::make_data_config(cmd);

    data_config.model_type = gelex::ModelType::A;
    data_config.grm_paths = assoc_config.grm_paths;
    data_config.transform_type = assoc_config.transform_type;
    data_config.int_offset = assoc_config.int_offset;

    int threads = cmd.get<int>("--threads");
    gelex::cli::setup_parallelization(threads);
    gelex::cli::print_assoc_header(threads);

    gelex::cli::DataPipeReporter data_reporter;
    gelex::DataPipe data(
        data_config,
        [&data_reporter](const gelex::DataPipeEvent& e)
        {
            std::visit([&](const auto& ev) { data_reporter.on_event(ev); }, e);
        });

    data.load();

    gelex::BedPipe bed_pipe(data_config.bed_path, data.sample_manager());
    auto bim = data_config.bed_path;
    auto snp_effects
        = std::move(gelex::detail::BimLoader(bim.replace_extension(".bim")))
              .take_info();

    gelex::cli::GwasRunner::Config runner_config{
        .max_iter = assoc_config.max_iter,
        .tol = assoc_config.tol,
        .chunk_size = data_config.chunk_size,
        .loco = assoc_config.loco,
        .additive = assoc_config.additive,
        .method = data_config.genotype_method,
        .grm_paths = assoc_config.grm_paths,
        .out_prefix = data_config.output_prefix};

    gelex::cli::GwasRunner runner(
        runner_config,
        std::move(data),
        std::move(bed_pipe),
        std::move(snp_effects));

    runner.run();

    return 0;
}
