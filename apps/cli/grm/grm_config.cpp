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

#include "grm_config.h"

#include <argparse.h>

#include "cli/cli_helper.h"
#include "gelex/data/genotype/bed_path.h"
#include "gelex/exception.h"
#include "gelex/types/freq_effect.h"

namespace gelex::cli
{

auto make_grm_config(argparse::ArgumentParser& cmd) -> gelex::GrmEngine::Config
{
    bool add = cmd.get<bool>("--add");
    bool dom = cmd.get<bool>("--dom");
    auto mode = gelex::freq::GrmType::A;
    if (add && dom)
    {
        mode = gelex::freq::GrmType::AD;
    }
    else if (dom)
    {
        mode = gelex::freq::GrmType::D;
    }
    auto chunk_size = cmd.get<int>("--chunk-size");
    if (chunk_size <= 0)
    {
        throw gelex::ArgumentValidationException("chunk_size must be positive");
    }

    return gelex::GrmEngine::Config{
        .bed_path = gelex::format_bed_path(cmd.get("--bfile")),
        .mode = mode,
        .method = gelex::cli::parse_genotype_process_method(
            cmd.get<std::string>("--geno-method")),
        .do_loco = cmd.get<bool>("--loco"),
        .out_prefix = cmd.get("--out"),
        .chunk_size = chunk_size};
}
}  // namespace gelex::cli
