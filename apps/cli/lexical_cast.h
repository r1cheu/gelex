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

#ifndef APPS_CLI_LEXICAL_CAST_H_
#define APPS_CLI_LEXICAL_CAST_H_

#include <string>

#include "gelex/data/genotype_method.h"
#include "gelex/infra/stats/rank_inverse_norm_transform.h"
#include "gelex/types/genetic_mode.h"

// Found by CLI11 via ADL, so add_option can bind enum targets directly.
namespace gelex
{
auto lexical_cast(const std::string& input, GenotypeMethod& output) -> bool;

auto lexical_cast(const std::string& input, GeneticMode& output) -> bool;

auto lexical_cast(const std::string& input, GeneticModeSet& output) -> bool;

auto lexical_cast(const std::string& input, RintType& output) -> bool;
}  // namespace gelex

#endif  // APPS_CLI_LEXICAL_CAST_H_
