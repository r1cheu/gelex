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

#include "cli/option_groups.h"

#include <CLI/CLI.hpp>

#include "cli/common_data.h"
#include "cli/random_design_data.h"
#include "cli/reml_data.h"
#include "cli/validators.h"

namespace cli
{

auto add_common_io_options(CLI::App& cmd, BaseDataConfig& config) -> void
{
    cmd.add_option(
           "-p,--pheno",
           config.pheno_path,
           "Phenotype TSV with FID, IID, trait columns")
        ->group("I/O")
        ->type_name("<PHENOTYPE>")
        ->check(CLI::ExistingFile)
        ->required();
    cmd.add_option(
           "--pheno-col",
           config.pheno_col,
           "0-based phenotype column after FID/IID; first trait=0")
        ->group("I/O")
        ->type_name("<COL>")
        ->check(cli::non_negative_number())
        ->capture_default_str();
    cmd.add_option(
           "--qcovar",
           config.qcovar_path,
           "Quantitative covariate TSV with FID, IID, numeric columns")
        ->group("I/O")
        ->type_name("<QCOVAR>")
        ->check(CLI::ExistingFile);
    cmd.add_option(
           "--dcovar",
           config.dcovar_path,
           "Discrete covariate TSV with FID, IID, factor columns")
        ->group("I/O")
        ->type_name("<DCOVAR>")
        ->check(CLI::ExistingFile);
}

auto add_random_design_options(CLI::App& cmd, RandomDesignDataConfig& config)
    -> void
{
    cmd.add_option(
           "--drand",
           config.drand_path,
           "Discrete random-effect TSV (FID\tIID\tFactor\t...); each factor "
           "column defines one one-hot random-effect block")
        ->group("I/O")
        ->type_name("<DRAND>")
        ->check(CLI::ExistingFile);
    cmd.add_option(
           "--qrand",
           config.qrand_paths,
           "Quantitative random-effect matrix TSV (FID\tIID\tValue\t...); "
           "each file defines one quantitative random-effect block")
        ->group("I/O")
        ->type_name("<QRAND>")
        ->expected(1, -1)
        ->allow_extra_args()
        ->check(CLI::ExistingFile);
}

auto add_random_effect_options(CLI::App& cmd, RemlDataConfig& config) -> void
{
    cmd.add_option(
           "--grm",
           config.grm,
           "GRM prefix(es) without suffix; reads <prefix>.bin/.id")
        ->group("I/O")
        ->type_name("<GRM>")
        ->expected(1, -1)
        ->allow_extra_args();
    add_random_design_options(cmd, config.designs);
    cmd.add_option(
           "--interaction",
           config.interactions,
           "Interaction random effect '<a>:<b>', the rescaled Hadamard product "
           "of two kernels. Each operand is a loaded effect name or a GRM "
           "prefix")
        ->group("I/O")
        ->type_name("<A:B>")
        ->expected(1, -1)
        ->allow_extra_args();
}

}  // namespace cli
