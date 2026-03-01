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

#include "assoc_config.h"

#include <argparse.h>

#include <filesystem>
#include <string_view>
#include <vector>

#include "gelex/data/genotype/bed_path.h"

namespace
{

auto parse_transform_type(std::string_view transform)
    -> gelex::detail::TransformType
{
    if (transform == "dint")
    {
        return gelex::detail::TransformType::DINT;
    }
    if (transform == "iint")
    {
        return gelex::detail::TransformType::IINT;
    }
    return gelex::detail::TransformType::None;
}

}  // namespace

auto AssocConfig::make(argparse::ArgumentParser& cmd) -> AssocConfig
{
    std::vector<std::filesystem::path> grm_paths;
    for (const auto& grm_path : cmd.get<std::vector<std::string>>("--grm"))
    {
        grm_paths.push_back(grm_path);
    }

    auto bed_path = gelex::format_bed_path(cmd.get("bfile"));
    auto transform_type = parse_transform_type(cmd.get("--transform"));

    auto genotype_method
        = gelex::parse_genotype_process_method(cmd.get<int>("--geno-method"));

    AssocConfig config{
        .threads = cmd.get<int>("--threads"),
        .phenotype_path = cmd.get("--pheno"),
        .phenotype_column = cmd.get<int>("--pheno-col"),
        .bed_path = bed_path,
        .chunk_size = cmd.get<int>("--chunk-size"),
        .qcovar_path = cmd.get("--qcovar"),
        .dcovar_path = cmd.get("--dcovar"),
        .out_prefix = cmd.get("--out"),
        .grm_paths = grm_paths,
        .transform_type = transform_type,
        .int_offset = cmd.get<double>("--int-offset"),
        .max_iter = cmd.get<int>("--max-iter"),
        .tol = cmd.get<double>("--tol"),
        .loco = cmd.get<bool>("--loco"),
        .additive = cmd.get("--model") == "a",
        .genotype_method = genotype_method};

    return config;
}
