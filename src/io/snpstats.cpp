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

#include "gelex/io/snpstats.h"

#include <fmt/format.h>
#include <Eigen/Core>
#include <cstdint>
#include <filesystem>
#include <utility>

#include "gelex/exception.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{

auto has_snp_stats(const BinaryReader& reader, GeneticMode mode) -> bool
{
    return reader.contains(fmt::format("{}/snp_stats", mode));
}

auto read_snp_stats(const BinaryReader& reader, GeneticMode mode) -> SnpStats
{
    auto stats_map = reader.to_map<double>(fmt::format("{}/snp_stats", mode));

    SnpStats data;
    data.code = stats_map.leftCols(3).transpose();
    data.A1freq = stats_map.col(3);

    if (const auto path = fmt::format("{}/geno_method", mode);
        reader.contains(path))
    {
        auto method_map = reader.to_map<uint8_t>(path);
        data.method = genotype_method_from_byte(method_map(0, 0));
    }

    if (const auto path = fmt::format("{}/valid_indices", mode);
        reader.contains(path))
    {
        auto valid_map = reader.to_map<int64_t>(path);
        const auto* src = valid_map.col(0).data();
        data.valid_indices.assign(src, src + valid_map.rows());
    }

    return data;
}

auto load_snp_stats(const std::filesystem::path& path) -> SnpStatsData
{
    BinaryReader reader(path.string());
    SnpStatsData data;
    for (const auto mode : ALL_GENETIC_MODES)
    {
        if (has_snp_stats(reader, mode))
        {
            data.emplace(mode, read_snp_stats(reader, mode));
        }
    }
    return data;
}

auto write_snp_stats(
    BinaryWriter& writer,
    GeneticMode mode,
    const SnpStats& stats) -> void
{
    const auto n_snps = stats.code.cols();
    if (stats.A1freq.size() != n_snps)
    {
        throw GelexException(
            fmt::format(
                "write_snp_stats: mismatched sizes code={}, A1freq={}",
                n_snps,
                stats.A1freq.size()));
    }

    writer.write(
        fmt::format("{}/geno_method", mode),
        static_cast<uint8_t>(std::to_underlying(stats.method)));

    const Eigen::MatrixXd code_t = stats.code.transpose();
    auto stats_handle
        = writer.reserve<double>(fmt::format("{}/snp_stats", mode), n_snps, 4);
    writer.write(stats_handle, code_t);
    writer.write(stats_handle, stats.A1freq);

    if (!stats.valid_indices.empty())
    {
        auto valid_handle = writer.reserve<int64_t>(
            fmt::format("{}/valid_indices", mode),
            stats.valid_indices.size(),
            1);
        writer.write(valid_handle, stats.valid_indices);
    }
}

}  // namespace gelex
