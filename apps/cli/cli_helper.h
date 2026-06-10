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

#ifndef GELEX_CLI_CLI_HELPER_H_
#define GELEX_CLI_CLI_HELPER_H_

#include <string>
#include <string_view>
#include <vector>

#include <argparse.h>
#include <barkeep.h>
#include <Eigen/Core>

#include "gelex/data/genotype/process_method.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::cli
{

auto parse_genotype_process_method(std::string_view value)
    -> GenotypeProcessMethod;

auto parse_genetic_modes(std::string_view sv) -> std::vector<GeneticMode>;

auto is_tty() -> bool;

auto setup_parallelization(int num_threads) -> void;

auto print_gelex_banner_message(std::string_view version) -> void;

auto format_epilog(std::string_view text) -> std::string;

}  // namespace gelex::cli

#endif  // GELEX_CLI_CLI_HELPER_H_
