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

#include "config.h"

#include <string>
#include <utility>
#include <vector>

#include <CLI/CLI.hpp>

#include "cli/cli_helper.h"
#include "gelex/engine/grm.h"
#include "gelex/exception.h"
#include "gelex/types/genetic_effect_type.h"

namespace cli
{

auto make_grm_config(CLI::App& cmd) -> gelex::GrmEngine::Config
{
    bool add = cmd.get_option("--add")->count() > 0;
    bool dom = cmd.get_option("--dom")->count() > 0;

    std::vector<gelex::GeneticMode> requested;
    if (add)
    {
        requested.push_back(gelex::GeneticMode::A);
    }
    if (dom)
    {
        requested.push_back(gelex::GeneticMode::D);
    }
    if (requested.empty())
    {
        requested.push_back(gelex::GeneticMode::A);
    }

    auto chunk_size = cmd.get_option("--chunk-size")->as<int>();
    if (chunk_size <= 0)
    {
        throw gelex::GelexException("chunk_size must be positive");
    }

    return gelex::GrmEngine::Config{
        .bfile_prefix = cmd.get_option("--bfile")->as<std::string>(),
        .requested_effects = std::move(requested),
        .method = cli::parse_genotype_method(
            cmd.get_option("--geno-method")->as<std::string>()),
        .do_loco = cmd.get_option("--loco")->count() > 0,
        .out_prefix = cmd.get_option("--out")->as<std::string>(),
        .chunk_size = chunk_size};
}
}  // namespace cli
