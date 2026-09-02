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

#include "gelex/data/dataframe/reader.h"

#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "gelex/data/reader.h"
#include "gelex/io/detail/parser.h"

namespace
{

auto detect_delimiter(const std::filesystem::path& path) -> char
{
    auto file = gelex::detail::open_file<std::ifstream>(path, std::ios::in);
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

auto snp_effects_schema(const std::filesystem::path& path)
    -> std::vector<gelex::ColumnType>
{
    std::ifstream file(path);
    std::string header_line;
    std::getline(file, header_line);

    std::vector<std::string> tokens;
    std::istringstream ss(header_line);
    std::string tok;
    while (std::getline(ss, tok, '\t'))
    {
        tokens.push_back(tok);
    }

    std::vector<gelex::ColumnType> schema;
    for (const auto& name : tokens)
    {
        if (name == "SNP")
        {
            continue;
        }
        auto type = (name == "CHR" || name == "A1" || name == "A2")
                        ? gelex::ColumnType::String
                        : gelex::ColumnType::Double;
        schema.push_back(type);
    }
    return schema;
}
}  // namespace

namespace gelex
{

auto read_fam(const std::filesystem::path& path) -> DataFrame<std::string>
{
    using enum ColumnType;
    constexpr std::array schema = {String, String, Int, String};
    ReadOptions options;
    options.header = false;
    options.delimiter = detect_delimiter(path);
    options.index_cols = {0, 1};
    options.names = {"father", "mother", "sex", "phenotype"};
    return read_dataframe<std::string>(path, options, schema);
}

auto read_bim(const std::filesystem::path& path) -> DataFrame<std::string>
{
    using enum ColumnType;
    constexpr std::array schema = {String, Int, String, String};
    ReadOptions options;
    options.header = false;
    options.delimiter = detect_delimiter(path);
    options.index_cols = {1};
    options.select_cols = {0, 3, 4, 5};
    options.names = {"CHR", "BP", "A1", "A2"};
    return read_dataframe<std::string>(path, options, schema);
}

auto read_snp_effects(const std::filesystem::path& path)
    -> DataFrame<std::string>
{
    ReadOptions options;
    options.index_cols = {1};
    return read_dataframe<std::string>(path, options, snp_effects_schema(path));
}

auto read_param(const std::filesystem::path& path) -> DataFrame<std::string>
{
    ReadOptions options;
    options.index_cols = {0};
    options.select_cols = {1};
    return read_dataframe<std::string, double>(path, options);
}

auto read_pheno(const std::filesystem::path& path, const std::size_t* pheno_col)
    -> DataFrame<std::string>
{
    ReadOptions options;
    options.index_cols = {0, 1};
    options.na_action = NaAction::Exclude;
    if (pheno_col != nullptr)
    {
        options.select_cols = {*pheno_col + 2};
    }
    return read_dataframe<std::string, double>(path, options);
}

auto read_qcovar(const std::filesystem::path& path) -> DataFrame<std::string>
{
    ReadOptions options;
    options.index_cols = {0, 1};
    return read_dataframe<std::string, double>(path, options);
}

auto read_dcovar(const std::filesystem::path& path) -> DataFrame<std::string>
{
    ReadOptions options;
    options.index_cols = {0, 1};
    return read_dataframe<std::string, std::string>(path, options);
}

}  // namespace gelex
