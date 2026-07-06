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

#include "gelex/io/snpstats/reader.h"

#include <fmt/format.h>
#include <cstdint>

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
    data.mean = stats_map.col(0);
    data.var = stats_map.col(1);
    data.A1freq = stats_map.col(2);

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

}  // namespace gelex
