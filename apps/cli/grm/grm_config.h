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

#ifndef GELEX_CLI_GRM_CONFIG_H_
#define GELEX_CLI_GRM_CONFIG_H_

#include <filesystem>
#include <string>

#include "gelex/data/genotype/genotype_processor.h"
#include "gelex/types/freq_effect.h"

namespace argparse
{
class ArgumentParser;
}

struct GrmConfig
{
    std::filesystem::path bed_path;
    std::string out_prefix;
    gelex::GenotypeProcessMethod method;
    gelex::freq::GrmType mode;
    int chunk_size;
    bool do_loco;
    int threads;

    static auto make(argparse::ArgumentParser& cmd) -> GrmConfig;

    auto validate() const -> void;
};

#endif  // GELEX_CLI_GRM_CONFIG_H_
