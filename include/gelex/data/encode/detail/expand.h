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

#ifndef GELEX_DATA_ENCODE_DETAIL_EXPAND_H_
#define GELEX_DATA_ENCODE_DETAIL_EXPAND_H_

#include <Eigen/Core>
#include <cstddef>
#include <cstdint>
#include <span>

namespace gelex::detail
{

// Expands one 2-bit-packed BED variant directly into its encoded column,
// gathering source samples to target rows through a 4-entry lookup indexed by
// the raw .bed code: 00 -> A1A1, 01 -> missing, 10 -> A1A2, 11 -> A2A2.
inline auto expand_encoded_column(
    std::span<const std::uint8_t> variant_bytes,
    std::span<const Eigen::Index> target_to_source,
    const Eigen::Ref<const Eigen::Array4d>& lut,
    Eigen::Ref<Eigen::VectorXd> out) -> void
{
    for (Eigen::Index row = 0;
         row < static_cast<Eigen::Index>(target_to_source.size());
         ++row)
    {
        const Eigen::Index source
            = target_to_source[static_cast<std::size_t>(row)];
        const auto byte = static_cast<unsigned>(
            variant_bytes[static_cast<std::size_t>(source / 4)]);
        const auto shift = static_cast<unsigned>(2 * (source % 4));
        out[row] = lut[(byte >> shift) & 0x03U];
    }
}

}  // namespace gelex::detail

#endif  // GELEX_DATA_ENCODE_DETAIL_EXPAND_H_
