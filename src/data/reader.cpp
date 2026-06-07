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
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <mio.h>
#include <Eigen/Core>

#include "gelex/data/dataframe/reader.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"
#include "gelex/io/detail/parser.h"

namespace
{

auto lower_triangle_index(size_t i, size_t j) -> size_t
{
    return (i * (i + 1) / 2) + j;
}

auto mmap_and_check_size(const std::string& path, size_t expected_size)
    -> mio::mmap_source
{
    std::error_code ec;
    mio::mmap_source mmap = mio::make_mmap_source(path, ec);
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
    const gelex::dataframe::Index<std::string>& source_indices,
    const gelex::dataframe::Index<std::string>& target_indices)
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
    auto file = gelex::io::detail::open_file<std::ifstream>(path, std::ios::in);
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
}  // namespace

namespace gelex
{

auto read_fam(const std::filesystem::path& path)
    -> dataframe::DataFrame<std::string>
{
    using enum dataframe::ColumnType;
    constexpr std::array SCHEMA = {String, String, Int, String};
    dataframe::ReadOptions options;
    options.header = false;
    options.delimiter = detect_delimiter(path);
    options.index_cols = {0, 1};
    options.names = {"father", "mother", "sex", "phenotype"};
    return dataframe::read_dataframe<std::string>(path, options, SCHEMA);
}

auto read_bim(const std::filesystem::path& path)
    -> dataframe::DataFrame<std::string>
{
    using enum dataframe::ColumnType;
    constexpr std::array SCHEMA = {String, Int, String, String};
    dataframe::ReadOptions options;
    options.header = false;
    options.delimiter = detect_delimiter(path);
    options.index_cols = {1};
    options.select_cols = {0, 3, 4, 5};
    options.names = {"chrom", "pos", "A1", "A2"};
    return dataframe::read_dataframe<std::string>(path, options, SCHEMA);
}

auto read_pheno(const std::filesystem::path& path, const std::size_t* pheno_col)
    -> dataframe::DataFrame<std::string>
{
    dataframe::ReadOptions options;
    options.index_cols = {0, 1};
    options.na_action = dataframe::NaAction::Exclude;
    if (pheno_col != nullptr)
    {
        options.select_cols = {*pheno_col + 2};
    }
    return dataframe::read_dataframe<std::string, double>(path, options);
}

auto read_qcovar(const std::filesystem::path& path)
    -> dataframe::DataFrame<std::string>
{
    dataframe::ReadOptions options;
    options.index_cols = {0, 1};
    return dataframe::read_dataframe<std::string, double>(path, options);
}

auto read_dcovar(const std::filesystem::path& path)
    -> dataframe::DataFrame<std::string>
{
    dataframe::ReadOptions options;
    options.index_cols = {0, 1};
    return dataframe::read_dataframe<std::string, std::string>(path, options);
}

auto read_grm_ids(const std::string& prefix)
    -> gelex::dataframe::Index<std::string>
{
    std::string path = prefix + ".id";
    dataframe::ReadOptions options;
    options.delimiter = '\t';
    options.header = false;
    options.index_cols = {0, 1};
    auto index = dataframe::read_index<std::string>(path, options);
    if (index.size() == 0)
    {
        throw GelexException(fmt::format("{}: no sample IDs found", path));
    }
    return index;
}

auto read_grm(
    const std::string& prefix,
    const dataframe::Index<std::string>* index,
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
