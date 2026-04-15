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

#include <array>
#include <fstream>
#include <string>

#include "gelex/data/dataframe/dataframe_reader.h"
#include "gelex/io/parser.h"

namespace gelex
{
inline auto detect_delimiter(const std::filesystem::path& path) -> char
{
    auto file = detail::open_file<std::ifstream>(path, std::ios::in);
    std::string line;
    std::getline(file, line);
    if (line.contains('\t'))
    {
        return '\t';
    }
    if (line.contains(','))
    {
        return ',';
    }
    if (line.contains(' '))
    {
        return ' ';
    }
    return '\t';
}

inline auto read_fam(const std::filesystem::path& path)
    -> df::DataFrame<std::string>
{
    using enum df::ColumnType;
    constexpr std::array kSchema = {String, String, Int, String};
    df::ReadOptions options;
    options.header = false;
    options.delimiter = detect_delimiter(path);
    options.index_cols = {0, 1};
    options.names = {"father", "mother", "sex", "phenotype"};
    return df::read_dataframe<std::string>(path, options, kSchema);
}

inline auto read_bim(const std::filesystem::path& path)
    -> df::DataFrame<std::string>
{
    using enum df::ColumnType;
    constexpr std::array kSchema = {String, Int, String, String};
    auto file = detail::open_file<std::ifstream>(path, std::ios::in);
    df::ReadOptions options;
    options.header = false;
    options.delimiter = detect_delimiter(path);
    options.index_cols = {1};
    options.select_cols = {0, 3, 4, 5};
    options.names = {"chrom", "pos", "A1", "A2"};
    return df::read_dataframe<std::string>(path, options, kSchema);
}

inline auto read_pheno(
    const std::filesystem::path& path,
    const std::size_t* pheno_col = nullptr) -> df::DataFrame<std::string>
{
    df::ReadOptions options;
    options.index_cols = {0, 1};
    options.na_action = df::NaAction::Exclude;
    if (pheno_col != nullptr)
    {
        options.select_cols = {*pheno_col + 2};
    }
    return df::read_dataframe<std::string, double>(path, options);
}

inline auto read_qcovar(const std::filesystem::path& path)
    -> df::DataFrame<std::string>
{
    df::ReadOptions options;
    options.index_cols = {0, 1};
    return df::read_dataframe<std::string, double>(path, options);
}

inline auto read_dcovar(const std::filesystem::path& path)
    -> df::DataFrame<std::string>
{
    df::ReadOptions options;
    options.index_cols = {0, 1};
    return df::read_dataframe<std::string, std::string>(path, options);
}

};  // namespace gelex

#endif  // GELEX_DATA_READER_H_
