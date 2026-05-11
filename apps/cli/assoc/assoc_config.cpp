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

#include "cli/cli_helper.h"
#include "gelex/algo/gwas/assoc_type.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::cli
{

auto parse_test_type(std::string_view test) -> AssocType
{
    if (test == "joint")
    {
        return AssocType::Joint;
    }
    return AssocType::Single;
}

auto parse_transform_type(std::string_view transform) -> detail::TransformType
{
    if (transform == "dint")
    {
        return detail::TransformType::DINT;
    }
    if (transform == "iint")
    {
        return detail::TransformType::IINT;
    }
    return detail::TransformType::None;
}

auto make_assoc_config(argparse::ArgumentParser& cmd) -> AssocEngine::Config
{
    auto test_type = parse_test_type(cmd.get("--test"));
    // Joint tester ignores mode; default to A for single when not specified
    GeneticMode mode = GeneticMode::A;
    if (test_type == AssocType::Single)
    {
        auto sv = cmd.get("--mode");
        if (sv == "D")
        {
            mode = GeneticMode::D;
        }
    }

    return AssocEngine::Config{
        .mode = mode,
        .method
        = parse_genotype_process_method(cmd.get<std::string>("--geno-method")),
        .chunk_size = cmd.get<int>("--chunk-size"),
        .max_iter = cmd.get<int>("--max-iter"),
        .tol = cmd.get<double>("--tol"),
        .bfile_prefix = cmd.get("--bfile"),
        .out_prefix = cmd.get("--out"),
        .loco = cmd.get<bool>("--loco"),
        .test_type = test_type};
}

}  // namespace gelex::cli
