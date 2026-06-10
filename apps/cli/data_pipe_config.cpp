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

#include <fmt/format.h>
#include <CLI/CLI.hpp>
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

namespace cli
{

auto make_pheno_config(CLI::App& cmd) -> gelex::PhenoPipe::Config
{
    return {
        .phenotype_path = cmd.get_option("--pheno")->as<std::string>(),
        .phenotype_column = cmd.get_option("--pheno-col")->as<int>(),
        .quantitative_covariates_path
        = cmd.get_option("--qcovar")->count() > 0
              ? std::make_optional(
                    std::filesystem::path(
                        cmd.get_option("--qcovar")->as<std::string>()))
              : std::nullopt,
        .discrete_covariates_path
        = cmd.get_option("--dcovar")->count() > 0
              ? std::make_optional(
                    std::filesystem::path(
                        cmd.get_option("--dcovar")->as<std::string>()))
              : std::nullopt,
    };
}

auto make_dataset_configs(CLI::App& cmd, bool use_mmap)
    -> std::pair<gelex::PhenoPipe::Config, gelex::GenoPipe::Config>
{
    gelex::PhenoPipe::Config pheno_config{
        .phenotype_path = cmd.get_option("--pheno")->as<std::string>(),
        .phenotype_column = cmd.get_option("--pheno-col")->as<int>(),
        .quantitative_covariates_path
        = cmd.get_option("--qcovar")->count() > 0
              ? std::make_optional(
                    std::filesystem::path(
                        cmd.get_option("--qcovar")->as<std::string>()))
              : std::nullopt,
        .discrete_covariates_path
        = cmd.get_option("--dcovar")->count() > 0
              ? std::make_optional(
                    std::filesystem::path(
                        cmd.get_option("--dcovar")->as<std::string>()))
              : std::nullopt,
    };

    gelex::GenoPipe::Config geno_config{
        .bfile_prefix = cmd.get_option("--bfile")->as<std::string>(),
        .requested_effects
        = parse_genetic_modes(cmd.get_option("--mode")->as<std::string>()),
        .genotype_method = parse_genotype_method(
            cmd.get_option("--geno-method")->as<std::string>()),
        .use_mmap = use_mmap,
        .chunk_size = cmd.get_option("--chunk-size")->as<int>(),
        .output_prefix = cmd.get_option("--out")->as<std::string>(),
    };

    return {std::move(pheno_config), std::move(geno_config)};
}

auto detail::intersect_or_throw_impl(
    std::vector<const gelex::dataframe::Index<std::string>*> indices,
    const gelex::DatasetObserver& observer,
    std::string_view what) -> gelex::dataframe::Index<std::string>
{
    auto common = gelex::dataframe::intersect<std::string>(indices);
    notify(observer, gelex::IntersectionEvent{.common_samples = common.size()});

    if (common.size() == 0)
    {
        throw gelex::GelexException(
            fmt::format(
                "No common samples across {}. Check that sample IDs match "
                "across input files.",
                what));
    }
    return common;
}

}  // namespace cli
