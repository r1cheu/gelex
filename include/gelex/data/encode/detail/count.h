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

#ifndef GELEX_DATA_ENCODE_DETAIL_COUNT_H_
#define GELEX_DATA_ENCODE_DETAIL_COUNT_H_

#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ranges>
#include <span>

#include "gelex/data/encode/detail/sample_mask.h"
#include "gelex/data/encode/stats.h"

namespace gelex::detail
{

// Counts the four genotype classes of one 2-bit-packed BED variant over the
// kept samples, fusing decode and tabulation via bit-plane popcounts. The mask
// carries both the sample subset and the sample/byte padding (its zero bits),
// so this is a single branchless pass regardless of subsetting. PLINK1 .bed
// spec:
//   00 -> nA1A1 (dosage 2)   01 -> n_missing
//   10 -> nA1A2 (dosage 1)   11 -> nA2A2 (dosage 0)
[[nodiscard]] inline auto count_genotypes(
    std::span<const std::uint8_t> variant_bytes,
    const SampleMask& mask) -> LocusStats
{
    static_assert(
        std::endian::native == std::endian::little,
        "count_genotypes assumes little-endian packed layout");

    constexpr std::uint64_t LO = SampleMask::SLOT_MASK;

    const auto keep_words = mask.words();
    const std::size_t num_bytes = variant_bytes.size();

    assert(
        num_bytes == static_cast<std::size_t>((mask.source_size() + 3) / 4)
        && "variant_bytes length must match mask sample count");

    LocusStats stats;

    for (const auto [word_idx, keep] : std::views::enumerate(keep_words))
    {
        const auto byte_offset = static_cast<std::size_t>(word_idx) * 8;
        const std::size_t take
            = std::min<std::size_t>(8, num_bytes - byte_offset);

        std::uint64_t genotype{0};
        std::memcpy(&genotype, variant_bytes.data() + byte_offset, take);

        const std::uint64_t lo = genotype & LO;
        const std::uint64_t hi = (genotype >> 1) & LO;

        stats.n_missing += std::popcount((lo & ~hi) & keep);  // 01
        stats.nA1A2 += std::popcount((~lo & hi) & keep);      // 10
        stats.nA2A2 += std::popcount((lo & hi) & keep);       // 11
    }

    stats.nA1A1 = mask.n_kept() - stats.n_missing - stats.nA1A2 - stats.nA2A2;

    return stats;
}

}  // namespace gelex::detail

#endif  // GELEX_DATA_ENCODE_DETAIL_COUNT_H_
