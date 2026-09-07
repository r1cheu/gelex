// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

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
