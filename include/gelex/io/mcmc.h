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

#ifndef GELEX_IO_MCMC_H_
#define GELEX_IO_MCMC_H_

#include <filesystem>
#include <string_view>

namespace gelex
{

class BayesModel;

namespace mcmc
{

class Result;

auto write_summary(const Result& result, std::string_view prefix) -> void;
auto write_params(const Result& result, std::string_view prefix) -> void;
auto write_snp_eff(
    const Result& result,
    const BayesModel& model,
    const std::filesystem::path& bim_path,
    std::string_view prefix) -> void;

}  // namespace mcmc

}  // namespace gelex

#endif  // GELEX_IO_MCMC_H_
