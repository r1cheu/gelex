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

#include "gelex/data/pipe/pheno.h"
#include "gelex/engine/assoc.h"

namespace argparse
{
class ArgumentParser;
}

namespace gelex::cli
{

auto make_assoc_config(argparse::ArgumentParser& cmd) -> AssocEngine::Config;

auto parse_transform_type(std::string_view transform)
    -> gelex::detail::TransformType;

}  // namespace gelex::cli

#endif  // GELEX_CLI_ASSOC_CONFIG_H_
