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

#include <filesystem>

#include <Eigen/Core>

#include "assoc_config.h"
#include "cli/cli_helper.h"
#include "cli/data_pipe_reporter.h"
#include "gelex/algo/infer/reml.h"
#include "gelex/data/genotype/bed_pipe.h"
#include "gelex/infra/logging/data_pipe_event.h"
#include "gelex/pipeline/data_pipe.h"
#include "gwas_runner.h"

#include "gelex/data/loader/bim_loader.h"

auto assoc_execute(argparse::ArgumentParser& cmd) -> int
{
    auto config = AssocConfig::make(cmd);

    gelex::cli::setup_parallelization(config.threads);
    gelex::cli::print_assoc_header(config.threads);

    gelex::DataPipe::Config data_pipe_config{
        .phenotype_path = config.phenotype_path,
        .phenotype_column = config.phenotype_column,
        .bed_path = config.bed_path,
        .use_dominance_effect = false,
        .use_mmap = false,
        .chunk_size = config.chunk_size,
        .qcovar_path = config.qcovar_path,
        .dcovar_path = config.dcovar_path,
        .output_prefix = config.out_prefix,
        .grm_paths = config.grm_paths,
        .transform_type = config.transform_type,
        .int_offset = config.int_offset};

    gelex::cli::DataPipeReporter pipe_reporter;
    auto data_pipe = gelex::load_data_for_reml(
        data_pipe_config,
        [&pipe_reporter](const gelex::DataPipeEvent& e)
        {
            std::visit([&](const auto& ev) { pipe_reporter.on_event(ev); }, e);
        });

    auto sample_manager = data_pipe.sample_manager();

    gelex::BedPipe bed_pipe(config.bed_path, sample_manager);
    auto bim_path = config.bed_path;
    bim_path.replace_extension(".bim");
    auto snp_effects
        = std::move(gelex::detail::BimLoader(bim_path)).take_info();

    gelex::cli::GwasRunner::Config runner_config{
        .max_iter = config.max_iter,
        .tol = config.tol,
        .chunk_size = config.chunk_size,
        .loco = config.loco,
        .additive = config.additive,
        .method = config.genotype_method,
        .grm_paths = config.grm_paths,
        .out_prefix = config.out_prefix};

    gelex::cli::GwasRunner runner(
        runner_config,
        std::move(data_pipe),
        std::move(bed_pipe),
        std::move(snp_effects));

    runner.run();

    return 0;
}
