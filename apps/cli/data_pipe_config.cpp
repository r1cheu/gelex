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

#include "data_pipe_config.h"

#include <argparse.h>
#include <optional>
#include <utility>

#include "cli/cli_helper.h"
#include "gelex/data/genotype/bed_path.h"

namespace gelex::cli
{

auto make_pheno_config(argparse::ArgumentParser& cmd) -> PhenoPipe::Config
{
    return {
        .phenotype_path = cmd.get("--pheno"),
        .phenotype_column = cmd.get<int>("--pheno-col"),
        .bed_path = format_bed_path(cmd.get<std::string>("--bfile")),
        .quantitative_covariates_path
        = cmd.is_used("--qcovar")
              ? std::make_optional(std::filesystem::path(cmd.get("--qcovar")))
              : std::nullopt,
        .discrete_covariates_path
        = cmd.is_used("--dcovar")
              ? std::make_optional(std::filesystem::path(cmd.get("--dcovar")))
              : std::nullopt,
    };
}

auto make_fit_data_configs(argparse::ArgumentParser& cmd, bool use_mmap)
    -> std::pair<PhenoPipe::Config, GenoPipe::Config>
{
    auto bed_path = format_bed_path(cmd.get<std::string>("--bfile"));

    PhenoPipe::Config pheno_config{
        .phenotype_path = cmd.get("--pheno"),
        .phenotype_column = cmd.get<int>("--pheno-col"),
        .bed_path = bed_path,
        .quantitative_covariates_path
        = cmd.is_used("--qcovar")
              ? std::make_optional(std::filesystem::path(cmd.get("--qcovar")))
              : std::nullopt,
        .discrete_covariates_path
        = cmd.is_used("--dcovar")
              ? std::make_optional(std::filesystem::path(cmd.get("--dcovar")))
              : std::nullopt,
    };

    GenoPipe::Config geno_config{
        .bed_path = bed_path,
        .genotype_method
        = parse_genotype_process_method(cmd.get<std::string>("--geno-method")),
        .use_mmap = use_mmap,
        .chunk_size = cmd.get<int>("--chunk-size"),
        .output_prefix = cmd.get("--out"),
    };

    return {std::move(pheno_config), std::move(geno_config)};
}

}  // namespace gelex::cli
