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

#include <algorithm>
#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ranges>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/data/dataframe/reader.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"
#include "gelex/io/detail/parser.h"
#include "gelex/io/mapped_file.h"

namespace
{

auto lower_triangle_index(size_t i, size_t j) -> size_t
{
    return (i * (i + 1) / 2) + j;
}

auto mmap_and_check_size(const std::string& path, size_t expected_size)
    -> gelex::MappedFile
{
    std::error_code ec;
    gelex::MappedFile mmap;
    mmap.map(path, ec);
    if (ec)
    {
        throw gelex::GelexException(
            fmt::format("failed to memory map {}.", path));
    }

    if (mmap.size() != expected_size)
    {
        throw gelex::GelexException(
            fmt::format(
                "{}: file size mismatch. Expected {} bytes, got {} bytes",
                path,
                expected_size,
                mmap.size()));
    }
    return mmap;
}

auto create_index_mapping(
    const gelex::DataFrameIndex<std::string>& source_indices,
    const gelex::DataFrameIndex<std::string>& target_indices)
    -> std::vector<std::pair<Eigen::Index, Eigen::Index>>
{
    std::vector<std::pair<Eigen::Index, Eigen::Index>> idx_mapping;
    idx_mapping.reserve(target_indices.size());

    for (auto&& [tgt_idx, id] : std::views::enumerate(target_indices.keys()))
    {
        if (!source_indices.contains(id))
        {
            throw gelex::GelexException(
                fmt::format(
                    "targe index '{}' not found in source indices", id));
        }

        idx_mapping.emplace_back(
            static_cast<Eigen::Index>(source_indices.at(id)),
            static_cast<Eigen::Index>(tgt_idx));
    }
    return idx_mapping;
}

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
    constexpr std::array SCHEMA = {String, String, Int, String};
    ReadOptions options;
    options.header = false;
    options.delimiter = detect_delimiter(path);
    options.index_cols = {0, 1};
    options.names = {"father", "mother", "sex", "phenotype"};
    return read_dataframe<std::string>(path, options, SCHEMA);
}

auto read_bim(const std::filesystem::path& path) -> DataFrame<std::string>
{
    using enum ColumnType;
    constexpr std::array SCHEMA = {String, Int, String, String};
    ReadOptions options;
    options.header = false;
    options.delimiter = detect_delimiter(path);
    options.index_cols = {1};
    options.select_cols = {0, 3, 4, 5};
    options.names = {"chrom", "pos", "A1", "A2"};
    return read_dataframe<std::string>(path, options, SCHEMA);
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

auto read_grm_ids(const std::string& prefix)
    -> gelex::DataFrameIndex<std::string>
{
    std::string path = prefix + ".id";
    ReadOptions options;
    options.delimiter = '\t';
    options.header = false;
    options.index_cols = {0, 1};
    auto index = read_index<std::string>(path, options);
    if (index.size() == 0)
    {
        throw GelexException(fmt::format("{}: no sample IDs found", path));
    }
    return index;
}

auto read_grm(
    const std::string& prefix,
    const DataFrameIndex<std::string>* index,
    bool normalize) -> Eigen::MatrixXd
{
    auto source_index = gelex::read_grm_ids(prefix);
    auto n = static_cast<Eigen::Index>(source_index.size());

    auto grm_path = prefix + ".bin";
    std::size_t buffer_size = n * (n + 1) / 2 * sizeof(float);
    auto mmap = mmap_and_check_size(grm_path, buffer_size);

    Eigen::MatrixXd target;

    if (index != nullptr)
    {
        auto out_n = static_cast<Eigen::Index>(index->size());
        auto idx_mapping = create_index_mapping(source_index, *index);
        target.setZero(out_n, out_n);

        const auto* data = reinterpret_cast<const float*>(mmap.data());

        for (Eigen::Index ii = 0; ii < out_n; ++ii)
        {
            auto [src_i, tgt_i] = idx_mapping[ii];

            for (Eigen::Index jj = 0; jj <= ii; ++jj)
            {
                auto [src_j, tgt_j] = idx_mapping[jj];

                // Read from lower triangle (ensure src_i >= src_j)
                Eigen::Index file_i = src_i;
                Eigen::Index file_j = src_j;
                if (file_i < file_j)
                {
                    std::swap(file_i, file_j);
                }

                auto idx = lower_triangle_index(file_i, file_j);
                target(tgt_i, tgt_j) = static_cast<double>(data[idx]);
                if (tgt_i != tgt_j)
                {
                    target(tgt_j, tgt_i) = target(tgt_i, tgt_j);
                }
            }
        }
    }
    else
    {
        target.resize(n, n);
        const auto* data = reinterpret_cast<const float*>(mmap.data());

        for (Eigen::Index i = 0; i < n; ++i)
        {
            for (Eigen::Index j = 0; j <= i; ++j)
            {
                auto idx = lower_triangle_index(i, j);
                target(i, j) = static_cast<double>(data[idx]);
                if (i != j)
                {
                    target(j, i) = target(i, j);  // Symmetric
                }
            }
        }
    }

    if (normalize && target.size() > 0)
    {
        double denominator
            = target.trace() / static_cast<double>(target.rows());
        target /= denominator;
    }

    return target;
}

}  // namespace gelex
