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

#ifndef GELEX_CLI_DATASET_LOADER_H_
#define GELEX_CLI_DATASET_LOADER_H_

#include <string>
#include <string_view>
#include <vector>

#include "gelex/data/dataframe/index.h"
#include "gelex/infra/logging/dataset_event.h"

namespace gelex::cli
{

// Compute sample intersection across the given indices, emit IntersectionEvent
// to `observer`, and throw GelexException if the result is empty. `what`
// describes the source set in the error message (e.g. "phenotype, genotype
// (.fam), and covariates").
auto intersect_or_throw(
    std::vector<const dataframe::Index<std::string>*> indices,
    const DatasetObserver& observer,
    std::string_view what) -> dataframe::Index<std::string>;

}  // namespace gelex::cli

#endif  // GELEX_CLI_DATASET_LOADER_H_
