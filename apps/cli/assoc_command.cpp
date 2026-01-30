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

#include <filesystem>
#include <memory>
#include <ranges>

#include <fmt/format.h>
#include <Eigen/Core>

#include "cli/cli_helper.h"
#include "cli/gwas_runner.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/data_pipe.h"
#include "gelex/estimator/freq/reml.h"
#include "gelex/logger.h"

#include "data/loader/bim_loader.h"

auto assoc_execute(argparse::ArgumentParser& cmd) -> int
{
    auto logger = gelex::logging::get();
    std::string out_prefix = cmd.get("--out");

    gelex::cli::setup_parallelization(cmd.get<int>("--threads"));

    auto grm_paths = std::ranges::to<std::vector<std::filesystem::path>>(
        cmd.get<std::vector<std::string>>("--grm"));

    auto bed_path = gelex::BedPipe::format_bed_path(cmd.get("bfile"));

    std::string transform_str = cmd.get("--transform");
    auto transform_type = gelex::detail::TransformType::None;
    if (transform_str == "dint")
    {
        transform_type = gelex::detail::TransformType::DINT;
    }
    else if (transform_str == "iint")
    {
        transform_type = gelex::detail::TransformType::IINT;
    }

    gelex::DataPipe::Config config{
        .phenotype_path = cmd.get("--pheno"),
        .phenotype_column = cmd.get<int>("--pheno-col"),
        .bed_path = bed_path,
        .use_dominance_effect = false,
        .use_mmap = false,
        .chunk_size = cmd.get<int>("--chunk-size"),
        .qcovar_path = cmd.get("--qcovar"),
        .dcovar_path = cmd.get("--dcovar"),
        .iid_only = cmd.get<bool>("--iid-only"),
        .output_prefix = cmd.get("--out"),
        .grm_paths = grm_paths,
        .transform_type = transform_type,
        .int_offset = cmd.get<double>("--int-offset")};

    gelex::cli::print_assoc_header(cmd.get<int>("--threads"));

    auto data_pipe = gelex::load_data_for_reml(config);

    auto sample_manager = data_pipe.sample_manager();

    gelex::BedPipe bed_pipe(bed_path, sample_manager);
    auto bim_path = bed_path;
    bim_path.replace_extension(".bim");
    auto snp_effects
        = std::move(gelex::detail::BimLoader(bim_path)).take_info();

    gelex::cli::GwasRunner::Config runner_config{
        .max_iter = cmd.get<int>("--max-iter"),
        .tol = cmd.get<double>("--tol"),
        .chunk_size = cmd.get<int>("--chunk-size"),
        .loco = cmd.get<bool>("--loco"),
        .additive = cmd.get("--model") == "a",
        .grm_paths = grm_paths,
        .out_prefix = out_prefix};

    gelex::cli::GwasRunner runner(
        runner_config,
        std::move(data_pipe),
        std::move(bed_pipe),
        std::move(snp_effects));

    runner.run();

    return 0;
}
