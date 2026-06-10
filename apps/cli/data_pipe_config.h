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

#ifndef APPS_CLI_DATA_PIPE_CONFIG_H_
#define APPS_CLI_DATA_PIPE_CONFIG_H_

#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "gelex/data/dataframe/index.h"
#include "gelex/data/pipe/geno.h"
#include "gelex/data/pipe/pheno.h"
#include "gelex/infra/logging/dataset_event.h"

namespace CLI
{
class App;
}

namespace cli
{

auto make_pheno_config(CLI::App& cmd) -> gelex::PhenoPipe::Config;

auto make_dataset_configs(CLI::App& cmd, bool use_mmap)
    -> std::pair<gelex::PhenoPipe::Config, gelex::GenoPipe::Config>;

namespace detail
{

auto intersect_or_throw_impl(
    std::vector<const gelex::dataframe::Index<std::string>*> indices,
    const gelex::DatasetObserver& observer,
    std::string_view what) -> gelex::dataframe::Index<std::string>;

}  // namespace detail

// Compute sample intersection across the given range sources, emit
// gelex::IntersectionEvent to `observer`, and throw gelex::GelexException if
// the result is empty. Each source must be a range whose element type is
// convertible to `const gelex::dataframe::Index<std::string>*`. `what`
// describes the source set in the error message (e.g. "phenotype, genotype
// (.fam), and covariates").
template <std::ranges::input_range... Sources>
    requires(
        std::convertible_to<
            std::ranges::range_reference_t<Sources>,
            const gelex::dataframe::Index<std::string>*>
        && ...)
auto intersect_or_throw(
    const gelex::DatasetObserver& observer,
    std::string_view what,
    Sources&&... sources) -> gelex::dataframe::Index<std::string>
{
    std::vector<const gelex::dataframe::Index<std::string>*> indices;
    (indices.append_range(std::forward<Sources>(sources)), ...);
    return detail::intersect_or_throw_impl(std::move(indices), observer, what);
}

}  // namespace cli

#endif  // APPS_CLI_DATA_PIPE_CONFIG_H_
