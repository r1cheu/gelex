// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_DATA_READER_H_
#define GELEX_DATA_READER_H_

#include <cstddef>
#include <filesystem>
#include <string>

#include "gelex/data/dataframe/dataframe.h"

namespace gelex
{
auto read_fam(const std::filesystem::path& path) -> DataFrame<std::string>;

auto read_bim(const std::filesystem::path& path) -> DataFrame<std::string>;

auto read_snp_effects(const std::filesystem::path& path)
    -> DataFrame<std::string>;

auto read_param(const std::filesystem::path& path) -> DataFrame<std::string>;

auto read_pheno(
    const std::filesystem::path& path,
    const std::size_t* pheno_col = nullptr) -> DataFrame<std::string>;

auto read_qcovar(const std::filesystem::path& path) -> DataFrame<std::string>;

auto read_dcovar(const std::filesystem::path& path) -> DataFrame<std::string>;

};  // namespace gelex

#endif  // GELEX_DATA_READER_H_
