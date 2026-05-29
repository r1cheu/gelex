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

#include "gelex/io/predict/input_reader.h"

#include <fmt/format.h>
#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <optional>
#include <span>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/dataframe/column.h"
#include "gelex/data/dataframe/constants.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/encode.h"
#include "gelex/data/dataframe/reader.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"

namespace gelex::predict
{

namespace detail
{

const std::unordered_set<std::string> kSnpStringCols{
    "Chrom",
    "Position",
    "A1",
    "A2"};

auto snp_effects_schema(const std::filesystem::path& path)
    -> std::vector<dataframe::ColumnType>
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

    std::vector<dataframe::ColumnType> schema;
    for (const auto& name : tokens)
    {
        if (name == "ID")
        {
            continue;
        }
        auto type = kSnpStringCols.contains(name)
                        ? dataframe::ColumnType::String
                        : dataframe::ColumnType::Double;
        schema.push_back(type);
    }
    return schema;
}

}  // namespace detail

auto read_coefficients(const std::filesystem::path& path) -> Coefficients
{
    dataframe::ReadOptions options;
    options.index_cols = {0};
    options.select_cols = {1};

    auto df = dataframe::read_dataframe<std::string, double>(path, options);
    return Coefficients{
        .names = std::move(df).index().take_keys(),
        .values = df["mean"].to_map<double>(),
    };
}

auto read_snp_effects(const std::filesystem::path& path)
    -> dataframe::DataFrame<std::string>
{
    dataframe::ReadOptions options;
    options.index_cols = {1};
    return dataframe::read_dataframe<std::string>(
        path, options, detail::snp_effects_schema(path));
}

auto read_covariates(
    const std::optional<std::filesystem::path>& qcovar_path,
    const std::optional<std::filesystem::path>& dcovar_path,
    const Coefficients& coefficients,
    dataframe::DataFrame<std::string>& sample_df) -> Eigen::MatrixXd
{
    // 1. read files
    std::optional<dataframe::DataFrame<std::string>> qcovar_df;
    std::optional<dataframe::DataFrame<std::string>> dcovar_df;

    if (qcovar_path)
    {
        qcovar_df = read_qcovar(*qcovar_path);
    }

    if (dcovar_path)
    {
        dcovar_df = read_dcovar(*dcovar_path);
    }

    // 2. intersect all DataFrames
    std::vector<dataframe::DataFrame<std::string>*> dfs;
    dfs.push_back(&sample_df);
    if (qcovar_df)
    {
        dfs.push_back(&*qcovar_df);
    }
    if (dcovar_df)
    {
        dfs.push_back(&*dcovar_df);
    }
    dataframe::intersect_inplace(
        std::span<dataframe::DataFrame<std::string>* const>{dfs});

    // 3. group dcovar terms by column name
    //    "Sex\x1FM" → col="Sex", level="M"
    std::map<std::string, std::vector<std::string>> dcovar_levels;
    for (const auto& name : coefficients.names)
    {
        auto pos = name.find(dataframe::kSeparator);
        if (pos != std::string::npos)
        {
            dcovar_levels[name.substr(0, pos)].push_back(name.substr(pos + 1));
        }
    }

    // 4. encode dcovar columns and check levels
    std::map<std::string, dataframe::EncodedResult<>> encoded;
    for (const auto& [col_name, levels] : dcovar_levels)
    {
        if (!dcovar_df)
        {
            throw GelexException(
                fmt::format(
                    "discrete covariate file required for term '{}'",
                    col_name));
        }
        auto& col = (*dcovar_df)[col_name];
        // TODO(rlchen): report level mismatch via observer
        [[maybe_unused]] auto mismatch = dataframe::check_levels(col, levels);
        encoded[col_name] = dataframe::encode(col, levels);
    }

    // 5. build design matrix in coefficients.names order
    const auto n_samples = static_cast<Eigen::Index>(sample_df.rows());
    const auto n_terms = static_cast<Eigen::Index>(coefficients.names.size());
    Eigen::MatrixXd X(n_samples, n_terms);

    for (Eigen::Index i = 0; i < n_terms; ++i)
    {
        const auto& term = coefficients.names[static_cast<std::size_t>(i)];

        if (term == "Intercept")
        {
            X.col(i).setOnes();
            continue;
        }

        auto sep_pos = term.find(dataframe::kSeparator);
        if (sep_pos != std::string::npos)
        {
            auto col_name = term.substr(0, sep_pos);
            const auto& enc = encoded.at(col_name);
            auto it = std::ranges::find(enc.level_names, term);
            auto col_idx = std::distance(enc.level_names.begin(), it);
            X.col(i) = enc.data.col(col_idx);
        }
        else
        {
            if (!qcovar_df)
            {
                throw GelexException(
                    fmt::format(
                        "quantitative covariate file required for term '{}'",
                        term));
            }
            X.col(i) = (*qcovar_df)[term].to_map<double>();
        }
    }
    return X;
}

}  // namespace gelex::predict
