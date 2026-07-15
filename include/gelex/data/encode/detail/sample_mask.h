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

#ifndef GELEX_DATA_ENCODE_DETAIL_SAMPLE_MASK_H_
#define GELEX_DATA_ENCODE_DETAIL_SAMPLE_MASK_H_

#include <Eigen/Core>
#include <cstdint>
#include <span>
#include <vector>

#include "gelex/data/detail/index_projection.h"

namespace gelex::detail
{

// Per-source-sample keep mask laid out to align with 2-bit-packed BED planes:
// source sample i occupies the low bit of its 2-bit slot at
// words_[i / 32], bit 2 * (i % 32). Odd bits are always zero, so a packed
// genotype plane can be masked with a single AND. Positions >= source_size
// (sample-count and byte padding) stay zero and are thereby excluded.
class SampleMask
{
   public:
    static constexpr std::uint64_t SLOT_MASK = 0x5555555555555555ULL;

    SampleMask(
        std::span<const Eigen::Index> target_to_source,
        Eigen::Index source_size);

    explicit SampleMask(const IndexProjection& projection)
        : SampleMask{projection.target_to_source(), projection.source_size()}
    {
    }

    [[nodiscard]] auto words() const noexcept -> std::span<const std::uint64_t>
    {
        return words_;
    }

    [[nodiscard]] auto source_size() const noexcept -> Eigen::Index
    {
        return source_size_;
    }

    [[nodiscard]] auto n_kept() const noexcept -> Eigen::Index
    {
        return n_kept_;
    }

   private:
    std::vector<std::uint64_t> words_;
    Eigen::Index source_size_{0};
    Eigen::Index n_kept_{0};
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_ENCODE_DETAIL_SAMPLE_MASK_H_
