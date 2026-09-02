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

#include "gelex/data/grm/io.h"

#include <Eigen/Core>
#include <algorithm>
#include <cstddef>
#include <fmt/format.h>
#include <ios>
#include <ranges>
#include <span>
#include <string>
#include <system_error>
#include <vector>

#include "gelex/data/dataframe/reader.h"
#include "gelex/data/sample_id.h"
#include "gelex/exception.h"
#include "gelex/io/detail/atomic_output_stream.h"
#include "gelex/io/detail/text_writer.h"
#include "gelex/io/mapped_file.h"

namespace
{

auto lower_triangle_index(std::size_t row, std::size_t column) -> std::size_t
{
    return (row * (row + 1) / 2) + column;
}

auto map_and_check_size(const std::string& path, std::size_t expected_size)
    -> gelex::MappedFile
{
    std::error_code error;
    gelex::MappedFile mapped_file;
    mapped_file.map(path, error);
    if (error)
    {
        throw gelex::GelexException(
            fmt::format("failed to memory map {}.", path));
    }

    if (mapped_file.size() != expected_size)
    {
        throw gelex::GelexException(
            fmt::format(
                "{}: file size mismatch. Expected {} bytes, got {} bytes",
                path,
                expected_size,
                mapped_file.size()));
    }
    return mapped_file;
}

struct IndexPosition
{
    Eigen::Index source;
    Eigen::Index target;
};

auto create_index_mapping(
    const gelex::DataFrameIndex<std::string>& source_index,
    std::span<const std::string> target_ids) -> std::vector<IndexPosition>
{
    std::vector<IndexPosition> mapping;
    mapping.reserve(target_ids.size());

    for (auto&& [target_position, id] : std::views::enumerate(target_ids))
    {
        if (!source_index.contains(id))
        {
            throw gelex::GelexException(
                fmt::format("target index '{}' not found in source index", id));
        }

        mapping.push_back(
            {.source = static_cast<Eigen::Index>(source_index.at(id)),
             .target = static_cast<Eigen::Index>(target_position)});
    }
    return mapping;
}

}  // namespace

namespace gelex
{
auto write_grm_ids(const std::string& prefix, std::span<const std::string> ids)
    -> void
{
    detail::TextWriter writer(prefix + ".id");
    for (const auto& id : ids)
    {
        auto [fid, iid] = split_sample_id(id);
        writer.write(fmt::format("{}\t{}", fid, iid));
    }
}

auto write_grm(
    const std::string& prefix,
    const Eigen::Ref<const Eigen::MatrixXd>& grm,
    std::span<const std::string> ids) -> void
{
    detail::AtomicOutputStream file(prefix + ".bin", std::ios::binary);

    auto n = grm.rows();
    auto m = grm.cols();

    if (n != m)
    {
        throw GelexException(
            fmt::format(
                "{}: GRM must be square, got {}x{}",
                file.path().string(),
                n,
                m));
    }
    if (n != static_cast<Eigen::Index>(ids.size()))
    {
        throw GelexException(
            fmt::format(
                "{}: Number of IDs ({}) does not match GRM size ({})",
                file.path().string(),
                ids.size(),
                n));
    }

    std::vector<float> row_buffer(static_cast<size_t>(n));
    for (Eigen::Index i = 0; i < n; ++i)
    {
        for (Eigen::Index j = 0; j <= i; ++j)
        {
            row_buffer[static_cast<size_t>(j)] = static_cast<float>(grm(i, j));
        }
        file.write(
            reinterpret_cast<const char*>(row_buffer.data()),
            static_cast<std::streamsize>((i + 1) * sizeof(float)));
    }
    file.commit();
    write_grm_ids(prefix, ids);
}

auto read_grm_ids(const std::string& prefix) -> DataFrameIndex<std::string>
{
    const std::string path = prefix + ".id";
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
    const auto source_index = read_grm_ids(prefix);
    const auto source_size = static_cast<Eigen::Index>(source_index.size());

    const std::string grm_path = prefix + ".bin";
    const auto element_count = static_cast<std::size_t>(source_size)
                               * static_cast<std::size_t>(source_size + 1) / 2;
    const auto mapped_file
        = map_and_check_size(grm_path, element_count * sizeof(float));
    const std::span<const float> data{
        reinterpret_cast<const float*>(mapped_file.data()), element_count};

    Eigen::MatrixXd target;
    if (index != nullptr)
    {
        const auto target_size = static_cast<Eigen::Index>(index->size());
        const auto mapping = create_index_mapping(source_index, index->keys());
        target.setZero(target_size, target_size);

        for (Eigen::Index target_row = 0; target_row < target_size;
             ++target_row)
        {
            const auto [source_row, output_row] = mapping[target_row];
            for (Eigen::Index target_column = 0; target_column <= target_row;
                 ++target_column)
            {
                const auto [source_column, output_column]
                    = mapping[target_column];
                const auto file_row = std::max(source_row, source_column);
                const auto file_column = std::min(source_row, source_column);
                const auto position = lower_triangle_index(
                    static_cast<std::size_t>(file_row),
                    static_cast<std::size_t>(file_column));
                target(output_row, output_column)
                    = static_cast<double>(data[position]);
                if (output_row != output_column)
                {
                    target(output_column, output_row)
                        = target(output_row, output_column);
                }
            }
        }
    }
    else
    {
        target.resize(source_size, source_size);
        for (Eigen::Index row = 0; row < source_size; ++row)
        {
            for (Eigen::Index column = 0; column <= row; ++column)
            {
                const auto position = lower_triangle_index(
                    static_cast<std::size_t>(row),
                    static_cast<std::size_t>(column));
                target(row, column) = static_cast<double>(data[position]);
                if (row != column)
                {
                    target.coeffRef(column, row) = target.coeff(row, column);
                }
            }
        }
    }

    if (normalize && target.size() > 0)
    {
        const double denominator
            = target.trace() / static_cast<double>(target.rows());
        target /= denominator;
    }
    return target;
}

}  // namespace gelex
