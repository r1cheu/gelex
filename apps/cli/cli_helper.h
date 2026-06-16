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

#ifndef APPS_CLI_CLI_HELPER_H_
#define APPS_CLI_CLI_HELPER_H_

#include <string_view>
#include <vector>

#include <barkeep.h>
#include <Eigen/Core>

#include "gelex/data/genotype_method.h"
#include "gelex/types/genetic_effect_type.h"

namespace CLI
{
class App;
}  // namespace CLI

namespace cli
{

auto parse_genotype_method(std::string_view value) -> gelex::GenotypeMethod;

auto parse_genetic_modes(std::string_view sv)
    -> std::vector<gelex::GeneticMode>;

auto is_tty() -> bool;

auto setup_parallelization(int num_threads) -> void;

auto add_common_io_options(CLI::App& cmd) -> void;

auto execute_cli_command(CLI::App& parser, int (*execute_fn)(CLI::App&)) -> int;

}  // namespace cli

#endif  // APPS_CLI_CLI_HELPER_H_
