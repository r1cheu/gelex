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
#include <fmt/format.h>
#include <filesystem>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "cli/cli_helper.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/pipe/geno.h"
#include "gelex/data/pipe/pheno.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/dataset_event.h"
#include "gelex/infra/logging/notify.h"

namespace gelex::cli
{

auto make_pheno_config(argparse::ArgumentParser& cmd) -> PhenoPipe::Config
{
    return {
        .phenotype_path = cmd.get("--pheno"),
        .phenotype_column = cmd.get<int>("--pheno-col"),
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

auto make_dataset_configs(argparse::ArgumentParser& cmd, bool use_mmap)
    -> std::pair<PhenoPipe::Config, GenoPipe::Config>
{
    PhenoPipe::Config pheno_config{
        .phenotype_path = cmd.get("--pheno"),
        .phenotype_column = cmd.get<int>("--pheno-col"),
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
        .bfile_prefix = cmd.get<std::string>("--bfile"),
        .requested_effects = parse_genetic_modes(cmd.get("--mode")),
        .genotype_method
        = parse_genotype_method(cmd.get<std::string>("--geno-method")),
        .use_mmap = use_mmap,
        .chunk_size = cmd.get<int>("--chunk-size"),
        .output_prefix = cmd.get("--out"),
    };

    return {std::move(pheno_config), std::move(geno_config)};
}

auto detail::intersect_or_throw_impl(
    std::vector<const dataframe::Index<std::string>*> indices,
    const DatasetObserver& observer,
    std::string_view what) -> dataframe::Index<std::string>
{
    auto common = dataframe::intersect<std::string>(indices);
    notify(observer, IntersectionEvent{.common_samples = common.size()});

    if (common.size() == 0)
    {
        throw GelexException(
            fmt::format(
                "No common samples across {}. Check that sample IDs match "
                "across input files.",
                what));
    }
    return common;
}

}  // namespace gelex::cli
