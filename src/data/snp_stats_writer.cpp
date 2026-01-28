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

#include "snp_stats_writer.h"

#include <array>
#include <cassert>
#include <cstdint>
#include <format>

#include "gelex/exception.h"
#include "parser.h"

namespace gelex::detail
{

SnpStatsWriter::SnpStatsWriter(const std::filesystem::path& file_path)
    : path_(file_path), io_buffer_(kDefaultBufferSize)
{
    file_ = detail::open_file<std::ofstream>(
        path_, std::ios::binary | std::ios::trunc, io_buffer_);
}

void SnpStatsWriter::write_data(
    const void* data,
    std::streamsize size,
    std::string_view error_msg)
{
    file_.write(reinterpret_cast<const char*>(data), size);
    if (!file_.good())
    {
        throw FileWriteException(
            std::format("{}:{}", path_.string(), error_msg));
    }
}

void SnpStatsWriter::write(
    int64_t num_samples,
    std::span<const int64_t> monomorphic_indices,
    std::span<const double> means,
    std::span<const double> stddevs)
{
    if (means.size() != stddevs.size())
    {
        throw ArgumentValidationException(
            std::format(
                "means ({}) and stddevs ({}) must have the same "
                "length.",
                means.size(),
                stddevs.size()));
    }

    if (means.empty() || stddevs.empty())
    {
        throw ArgumentValidationException("means and stddevs cannot be empty");
    }

    const auto num_variants = static_cast<int64_t>(means.size());
    const auto num_monomorphic
        = static_cast<int64_t>(monomorphic_indices.size());

    const std::array<int64_t, 3> header
        = {num_samples, num_variants, num_monomorphic};
    check_monomorphic_indices(monomorphic_indices, num_variants);
    write_data(std::span<const int64_t>(header), "failed to write header");

    if (!monomorphic_indices.empty())
    {
        write_data(
            monomorphic_indices, "failed to write monomorphic snp indices");
    }

    write_data(means, "failed to write means");

    write_data(stddevs, "failed to write stddevs");
}

void SnpStatsWriter::check_monomorphic_indices(
    std::span<const int64_t> monomorphic_indices,
    int64_t num_variants)
{
    if (monomorphic_indices.empty())
    {
        return;
    }
    const int64_t max_value = monomorphic_indices.back();
    if (max_value >= num_variants || max_value < 0)
    {
        throw ArgumentValidationException(
            std::format(
                "Monomorphic SNP index {} is out of range [0, {}).",
                max_value,
                num_variants));
    }
}

}  // namespace gelex::detail
