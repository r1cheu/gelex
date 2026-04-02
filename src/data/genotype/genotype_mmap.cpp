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

#include "gelex/data/genotype/genotype_mmap.h"

#include <fmt/format.h>
#include <algorithm>
#include <cstdint>
#include <filesystem>

#include "gelex/io/binary_reader.h"

namespace gelex
{

GenotypeMap::GenotypeMap(
    const std::filesystem::path& bin_file,
    GeneticKind effect_type)
    : mat_(nullptr, 0, 0)
{
    const auto effect = EffectType::from_genetic(effect_type);
    reader_ = std::make_unique<detail::BinaryReader>(bin_file.string());

    auto geno_map = reader_->to_map<double>(fmt::format("{}/genotype", effect));
    rows_ = geno_map.rows();
    cols_ = geno_map.cols();
    new (&mat_) MapType(geno_map.data(), rows_, cols_);

    auto stats_mat
        = reader_->to_mat<double>(fmt::format("{}/loci_stats", effect));
    mean_ = stats_mat.col(0);
    stddev_ = stats_mat.col(1);

    if (const auto path = fmt::format("{}/mono_indices", effect);
        reader_->contains(path))
    {
        auto mono_mat = reader_->to_mat<int64_t>(path);
        const auto* src = mono_mat.col(0).data();
        mono_indices_.assign(src, src + mono_mat.rows());
        std::ranges::sort(mono_indices_);
    }
}

bool GenotypeMap::is_monomorphic(Eigen::Index snp_index) const noexcept
{
    return std::ranges::binary_search(mono_indices_, snp_index);
}

}  // namespace gelex
