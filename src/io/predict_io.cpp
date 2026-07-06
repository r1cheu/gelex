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

#include "gelex/io/predict_io.h"

#include <fmt/format.h>
#include <Eigen/Core>
#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <optional>
#include <ranges>
#include <span>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "gelex/data/dataframe/column.h"
#include "gelex/data/dataframe/constants.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/encode.h"
#include "gelex/data/dataframe/reader.h"
#include "gelex/data/reader.h"
#include "gelex/data/sample_id.h"
#include "gelex/exception.h"
#include "gelex/io/detail/text_writer.h"
#include "gelex/predict/types.h"

namespace gelex
{

namespace detail
{

auto snp_effects_schema(const std::filesystem::path& path)
    -> std::vector<ColumnType>
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

    std::vector<ColumnType> schema;
    for (const auto& name : tokens)
    {
        if (name == "SNP")
        {
            continue;
        }
        auto type = (name == "CHR" || name == "A1" || name == "A2")
                        ? ColumnType::String
                        : ColumnType::Double;
        schema.push_back(type);
    }
    return schema;
}

}  // namespace detail

auto read_coefficients(const std::filesystem::path& path) -> Coefficients
{
    ReadOptions options;
    options.index_cols = {0};
    options.select_cols = {1};

    auto df = read_dataframe<std::string, double>(path, options);
    return Coefficients{
        .names = std::move(df).index().take_keys(),
        .values = df["mean"].to_map<double>(),
    };
}

auto read_snp_effects(const std::filesystem::path& path)
    -> DataFrame<std::string>
{
    ReadOptions options;
    options.index_cols = {1};
    return read_dataframe<std::string>(
        path, options, detail::snp_effects_schema(path));
}

auto read_covariates(
    const std::optional<std::filesystem::path>& qcovar_path,
    const std::optional<std::filesystem::path>& dcovar_path,
    const Coefficients& coefficients,
    DataFrame<std::string>& sample_df) -> Eigen::MatrixXd
{
    // 1. read files
    std::optional<DataFrame<std::string>> qcovar_df;
    std::optional<DataFrame<std::string>> dcovar_df;

    if (qcovar_path)
    {
        qcovar_df = read_qcovar(*qcovar_path);
    }

    if (dcovar_path)
    {
        dcovar_df = read_dcovar(*dcovar_path);
    }

    // 2. intersect all DataFrames
    std::vector<DataFrame<std::string>*> dfs;
    dfs.push_back(&sample_df);
    if (qcovar_df)
    {
        dfs.push_back(&*qcovar_df);
    }
    if (dcovar_df)
    {
        dfs.push_back(&*dcovar_df);
    }
    intersect_inplace(std::span<DataFrame<std::string>* const>{dfs});

    // 3. group dcovar terms by column name
    // "Sex\x1FM" → col="Sex", level="M"
    std::map<std::string, std::vector<std::string>> dcovar_levels;
    for (const auto& name : coefficients.names)
    {
        auto pos = name.find(SEPARATOR);
        if (pos != std::string::npos)
        {
            dcovar_levels[name.substr(0, pos)].push_back(name.substr(pos + 1));
        }
    }

    // 4. encode dcovar columns and check levels
    std::map<std::string, EncodedResult<>> encoded;
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
        [[maybe_unused]] auto mismatch = check_levels(col, levels);
        encoded[col_name] = encode(col, levels);
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

        auto sep_pos = term.find(SEPARATOR);
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

auto write_predictions(
    const std::filesystem::path& output_path,
    const PredictResult& result) -> void
{
    if (output_path.empty())
    {
        throw GelexException("Output path must be provided");
    }

    const auto n_samples = static_cast<Eigen::Index>(result.sample_ids.size());
    if (n_samples != result.predictions.size())
    {
        throw GelexException(
            fmt::format(
                "Dimension mismatch: {} sample IDs but {} predictions",
                result.sample_ids.size(),
                result.predictions.size()));
    }

    detail::TextWriter writer(output_path);

    std::string header = "FID\tIID\tprediction";
    for (const auto& name : result.covar_names)
    {
        header += '\t';
        header += name;
    }
    for (const auto& mode : std::views::keys(result.snp_components))
    {
        header += fmt::format("\t{}", mode);
    }
    writer.write(header);

    std::string row_buf;
    for (Eigen::Index i = 0; i < n_samples; ++i)
    {
        row_buf.clear();
        auto [fid, iid]
            = split_sample_id(result.sample_ids[static_cast<std::size_t>(i)]);
        row_buf += fmt::format("{}\t{}", fid, iid);
        row_buf += fmt::format("\t{:.6f}", result.predictions[i]);

        for (Eigen::Index j = 0; j < result.covar_predictions.cols(); ++j)
        {
            row_buf += fmt::format("\t{:.6f}", result.covar_predictions(i, j));
        }
        for (const auto& component : std::views::values(result.snp_components))
        {
            row_buf += fmt::format("\t{:.6f}", component[i]);
        }
        writer.write(row_buf);
    }
}

}  // namespace gelex
