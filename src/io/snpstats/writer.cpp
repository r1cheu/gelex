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

#include "gelex/io/snpstats/writer.h"

#include <fmt/format.h>
#include <Eigen/Core>
#include <cstdint>
#include <utility>

#include "gelex/exception.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{

auto write_snp_stats(
    BinaryWriter& writer,
    GeneticMode mode,
    const SnpStats& stats) -> void
{
    const auto n_snps = stats.mean.size();
    if (stats.var.size() != n_snps || stats.A1freq.size() != n_snps)
    {
        throw GelexException(
            fmt::format(
                "write_snp_stats: mismatched sizes mean={}, var={}, A1freq={}",
                n_snps,
                stats.var.size(),
                stats.A1freq.size()));
    }

    writer.write(
        fmt::format("{}/geno_method", mode),
        static_cast<uint8_t>(std::to_underlying(stats.method)));

    auto stats_handle
        = writer.reserve<double>(fmt::format("{}/snp_stats", mode), n_snps, 3);
    writer.write(stats_handle, stats.mean);
    writer.write(stats_handle, stats.var);
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
