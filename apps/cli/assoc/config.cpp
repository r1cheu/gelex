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

#include <CLI/CLI.hpp>
#include <string>
#include <string_view>

#include "cli/cli_helper.h"
#include "gelex/algo/gwas/assoc_type.h"
#include "gelex/data/pipe/pheno.h"
#include "gelex/engine/assoc.h"
#include "gelex/types/genetic_effect_type.h"

namespace cli
{

auto parse_test_type(std::string_view test) -> gelex::AssocType
{
    if (test == "joint")
    {
        return gelex::AssocType::Joint;
    }
    return gelex::AssocType::Single;
}

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

auto make_assoc_config(CLI::App& cmd) -> gelex::AssocEngine::Config
{
    auto test_type
        = parse_test_type(cmd.get_option("--test")->as<std::string>());
    // Joint tester ignores mode; default to A for single when not specified
    gelex::GeneticMode mode = gelex::GeneticMode::A;
    if (test_type == gelex::AssocType::Single)
    {
        auto sv = cmd.get_option("--mode")->as<std::string>();
        if (sv == "D")
        {
            mode = gelex::GeneticMode::D;
        }
    }

    return gelex::AssocEngine::Config{
        .mode = mode,
        .method = parse_genotype_method(
            cmd.get_option("--geno-method")->as<std::string>()),
        .chunk_size = cmd.get_option("--chunk-size")->as<int>(),
        .max_iter = cmd.get_option("--max-iter")->as<int>(),
        .tol = cmd.get_option("--tol")->as<double>(),
        .bfile_prefix = cmd.get_option("--bfile")->as<std::string>(),
        .out_prefix = cmd.get_option("--out")->as<std::string>(),
        .loco = cmd.get_option("--loco")->count() > 0,
        .test_type = test_type};
}

}  // namespace cli
