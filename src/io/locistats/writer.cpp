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

#include "gelex/io/locistats/writer.h"

#include <fmt/format.h>
#include <Eigen/Core>
#include <cstdint>
#include <span>
#include <string_view>

#include "gelex/exception.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{

LociStatsWriter::LociStatsWriter(std::string_view output_path)
    : writer_(output_path)
{
}

auto LociStatsWriter::write(
    GeneticMode mode,
    uint8_t method,
    const Eigen::VectorXd& mean,
    const Eigen::VectorXd* stddev,
    std::span<const int64_t> mono_indices) -> void
{
    if (stddev != nullptr && mean.size() != stddev->size())
    {
        throw GelexException(
            fmt::format(
                "LociStatsWriter: mean size ({}) != stddev size ({})",
                mean.size(),
                stddev->size()));
    }

    const auto n_snps = mean.size();
    const Eigen::Index n_cols = (stddev != nullptr) ? 2 : 1;

    auto handle
        = writer_.reserve<uint8_t>(fmt::format("{}/geno_method", mode), 1, 1);
    writer_.write(handle, method);

    auto stats_handle = writer_.reserve<double>(
        fmt::format("{}/loci_stats", mode), n_snps, n_cols);

    writer_.write(stats_handle, mean);

    if (stddev != nullptr)
    {
        writer_.write(stats_handle, *stddev);
    }

    if (!mono_indices.empty())
    {
        auto mono_handle = writer_.reserve<int64_t>(
            fmt::format("{}/mono_indices", mode), mono_indices.size(), 1);
        writer_.write(mono_handle, mono_indices);
    }
}

}  // namespace gelex
