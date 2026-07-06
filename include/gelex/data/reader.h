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

#ifndef GELEX_DATA_READER_H_
#define GELEX_DATA_READER_H_

#include <cstddef>
#include <filesystem>
#include <string>

#include <Eigen/Core>

#include "gelex/data/dataframe/dataframe.h"

namespace gelex
{
auto read_fam(const std::filesystem::path& path) -> DataFrame<std::string>;

auto read_bim(const std::filesystem::path& path) -> DataFrame<std::string>;

auto read_pheno(
    const std::filesystem::path& path,
    const std::size_t* pheno_col = nullptr) -> DataFrame<std::string>;

auto read_qcovar(const std::filesystem::path& path) -> DataFrame<std::string>;

auto read_dcovar(const std::filesystem::path& path) -> DataFrame<std::string>;

auto read_grm_ids(const std::string& prefix)
    -> gelex::DataFrameIndex<std::string>;

auto read_grm(
    const std::string& prefix,
    const DataFrameIndex<std::string>* index = nullptr,
    bool normalize = true) -> Eigen::MatrixXd;

};  // namespace gelex

#endif  // GELEX_DATA_READER_H_
