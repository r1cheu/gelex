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

auto GrmConfig::make(argparse::ArgumentParser& cmd) -> GrmConfig
{
    bool add = cmd.get<bool>("--add");
    bool dom = cmd.get<bool>("--dom");
    gelex::freq::GrmType mode = gelex::freq::GrmType::A;
    if (add && dom)
    {
        mode = gelex::freq::GrmType::AD;
    }
    else if (dom)
    {
        mode = gelex::freq::GrmType::D;
    }

    return GrmConfig{
        .bed_path = gelex::format_bed_path(cmd.get("--bfile")),
        .out_prefix = cmd.get("--out"),
        .method = gelex::cli::parse_genotype_process_method(
            cmd.get<std::string>("--geno-method")),
        .mode = mode,
        .chunk_size = cmd.get<int>("--chunk-size"),
        .do_loco = cmd.get<bool>("--loco"),
        .threads = cmd.get<int>("--threads")};
}

auto GrmConfig::validate() const -> void
{
    if (chunk_size <= 0)
    {
        throw gelex::ArgumentValidationException("chunk_size must be positive");
    }
    if (threads == 0 || threads < -1)
    {
        throw gelex::ArgumentValidationException(
            "threads must be -1 (auto) or a positive integer");
    }
}
