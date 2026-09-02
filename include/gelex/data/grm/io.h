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

#ifndef GELEX_DATA_GRM_IO_H_
#define GELEX_DATA_GRM_IO_H_

#include <Eigen/Dense>
#include <span>
#include <string>

#include "gelex/data/dataframe/index.h"

namespace gelex
{
auto write_grm_ids(const std::string& prefix, std::span<const std::string> ids)
    -> void;

auto write_grm(
    const std::string& prefix,
    const Eigen::Ref<const Eigen::MatrixXd>& grm,
    std::span<const std::string> ids) -> void;

auto read_grm_ids(const std::string& prefix) -> DataFrameIndex<std::string>;

auto read_grm(
    const std::string& prefix,
    const DataFrameIndex<std::string>* index = nullptr,
    bool normalize = true) -> Eigen::MatrixXd;

}  // namespace gelex

#endif  // GELEX_DATA_GRM_IO_H_
