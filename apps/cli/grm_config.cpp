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

#include "gelex/data/genotype/bed_path.h"

auto GrmConfig::make(argparse::ArgumentParser& cmd) -> GrmConfig
{
    GrmConfig config{
        .bed_path = gelex::format_bed_path(cmd.get("--bfile")),
        .out_prefix = cmd.get("--out"),
        .method = cmd.get("--geno-method"),
        .chunk_size = cmd.get<int>("--chunk-size"),
        .do_additive = cmd.get<bool>("--add"),
        .do_dominant = cmd.get<bool>("--dom"),
        .do_loco = cmd.get<bool>("--loco"),
        .threads = cmd.get<int>("--threads")};

    if (!config.do_additive && !config.do_dominant)
    {
        config.do_additive = true;
    }

    return config;
}
