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

#include <string_view>
#include <vector>

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

namespace gelex::cli
{

auto make_assoc_config(argparse::ArgumentParser& cmd) -> AssocConfig
{
    std::vector<std::filesystem::path> grm_paths;
    for (const auto& p : cmd.get<std::vector<std::string>>("--grm"))
    {
        grm_paths.push_back(p);
    }

    return AssocConfig{
        .max_iter = cmd.get<int>("--max-iter"),
        .tol = cmd.get<double>("--tol"),
        .loco = cmd.get<bool>("--loco"),
        .additive = cmd.get("--model") == "a",
        .grm_paths = std::move(grm_paths),
        .transform_type = parse_transform_type(cmd.get("--transform")),
        .int_offset = cmd.get<double>("--int-offset")};
}

}  // namespace gelex::cli
