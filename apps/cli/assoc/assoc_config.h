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

#ifndef GELEX_CLI_ASSOC_CONFIG_H_
#define GELEX_CLI_ASSOC_CONFIG_H_

#include <filesystem>
#include <string>
#include <vector>

#include "gelex/data/genotype/genotype_method_dispatch.h"
#include "gelex/pipeline/data_pipe.h"

namespace argparse
{
class ArgumentParser;
}

struct AssocConfig
{
    int threads;
    std::filesystem::path phenotype_path;
    int phenotype_column;
    std::filesystem::path bed_path;
    int chunk_size;
    std::filesystem::path qcovar_path;
    std::filesystem::path dcovar_path;
    std::string out_prefix;
    std::vector<std::filesystem::path> grm_paths;
    gelex::detail::TransformType transform_type;
    double int_offset;
    int max_iter;
    double tol;
    bool loco;
    bool additive;
    gelex::GenotypeProcessMethod genotype_method;

    static auto make(argparse::ArgumentParser& cmd) -> AssocConfig;
};

#endif  // GELEX_CLI_ASSOC_CONFIG_H_
