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

#include "gelex/data/grm/grm_reader.h"

#include <fmt/format.h>
#include <ranges>
#include <string>
#include <system_error>

#include "gelex/data/dataframe/dataframe_reader.h"
#include "gelex/exception.h"
#include "gelex/types/freq_effect.h"

namespace
{
auto get_type(std::string_view grm_path_stem) -> gelex::GeneticMode
{
    if (grm_path_stem.contains("add"))
    {
        return gelex::GeneticMode::A;
    }
    if (grm_path_stem.contains("dom"))
    {
        return gelex::GeneticMode::D;
    }
    return gelex::GeneticMode::A;
}
}  // namespace
namespace gelex::detail
{

namespace
{
auto read_grm_index(const std::filesystem::path& id_path)
    -> gelex::df::Index<std::string>
{
    df::ReadOptions options;
    options.delimiter = '\t';
    options.header = false;
    options.index_cols = {0, 1};
    auto index = df::read_index<std::string>(id_path, options);
    if (index.size() == 0)
    {
        throw GelexException(
            fmt::format("{}: no sample IDs found", id_path.string()));
    }
    return index;
}
}  // namespace

GrmReader::GrmReader(const std::filesystem::path& prefix)
    : bin_path_(prefix.string() + ".bin"),
      sample_index_(read_grm_index(prefix.string() + ".id")),
      type_(get_type(prefix.string()))
{
    init_mmap();
}

auto GrmReader::init_mmap() -> void
{
    std::error_code ec;
    mmap_.map(bin_path_.string(), ec);
    if (ec)
    {
        throw GelexException(
            fmt::format("{}: failed to mmap file", bin_path_.string()));
    }

    // GRM binary format: [float32 lower triangle]
    // Expected size = n * (n + 1) / 2 * sizeof(float)
    auto n = static_cast<size_t>(num_samples());
    size_t expected_elements = n * (n + 1) / 2;
    size_t expected_size = expected_elements * sizeof(float);

    if (mmap_.size() != expected_size)
    {
        throw GelexException(
            fmt::format(
                "{}: file size mismatch. Expected {} bytes ({} samples), got "
                "{} bytes",
                bin_path_.string(),
                expected_size,
                num_samples(),
                mmap_.size()));
    }
}

auto GrmReader::load() const -> Eigen::MatrixXd
{
    Eigen::MatrixXd grm = load_unnormalized();
    double denominator = grm.trace() / static_cast<double>(num_samples());
    grm /= denominator;
    return grm;
}

auto GrmReader::load_unnormalized() const -> Eigen::MatrixXd
{
    auto n = num_samples();
    Eigen::MatrixXd grm(n, n);

    const auto* data = reinterpret_cast<const float*>(mmap_.data());

    // Fill lower triangle and mirror to upper triangle
    for (Eigen::Index i = 0; i < n; ++i)
    {
        for (Eigen::Index j = 0; j <= i; ++j)
        {
            size_t idx = lower_triangle_index(i, j);
            auto value = static_cast<double>(data[idx]);
            grm(i, j) = value;
            grm(j, i) = value;
        }
    }

    return grm;
}

auto GrmReader::load(const df::Index<std::string>& sample_index) const
    -> Eigen::MatrixXd
{
    Eigen::MatrixXd grm = load_unnormalized(sample_index);
    if (grm.size() > 0)
    {
        double denominator = grm.trace() / static_cast<double>(grm.rows());
        grm /= denominator;
    }
    return grm;
}

auto GrmReader::load_unnormalized(
    const df::Index<std::string>& sample_index,
    Eigen::MatrixXd& target) const -> void
{
    if (sample_index.size() == 0)
    {
        target.resize(0, 0);
        return;
    }

    std::vector<std::pair<Eigen::Index, Eigen::Index>> idx_mapping;
    idx_mapping.reserve(sample_index.size());

    for (auto&& [tgt_idx, id] : std::views::enumerate(sample_index.keys()))
    {
        if (!sample_index_.contains(id))
        {
            throw GelexException(
                fmt::format(
                    "{}: sample ID '{}' not found in GRM file",
                    bin_path_.string(),
                    id));
        }

        idx_mapping.emplace_back(
            static_cast<Eigen::Index>(sample_index_.at(id)),
            static_cast<Eigen::Index>(tgt_idx));
    }

    // Allocate output matrix
    auto out_size = static_cast<Eigen::Index>(sample_index.size());
    target.setZero(out_size, out_size);

    const auto* data = reinterpret_cast<const float*>(mmap_.data());

    // Fill the matrix using index mapping
    for (size_t ii = 0; ii < idx_mapping.size(); ++ii)
    {
        auto [src_i, tgt_i] = idx_mapping[ii];

        for (size_t jj = 0; jj < idx_mapping.size(); ++jj)
        {
            auto [src_j, tgt_j] = idx_mapping[jj];

            // Read from lower triangle (ensure src_i >= src_j)
            Eigen::Index file_i = src_i;
            Eigen::Index file_j = src_j;
            if (file_i < file_j)
            {
                std::swap(file_i, file_j);
            }

            size_t idx = lower_triangle_index(file_i, file_j);
            target(tgt_i, tgt_j) = static_cast<double>(data[idx]);
        }
    }
}

auto GrmReader::load_unnormalized(
    const df::Index<std::string>& sample_index) const -> Eigen::MatrixXd
{
    Eigen::MatrixXd grm;
    load_unnormalized(sample_index, grm);
    return grm;
}

}  // namespace gelex::detail
