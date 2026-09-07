// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/data/encode/detail/sample_mask.h"

#include <cstddef>

namespace gelex::detail
{

namespace
{
constexpr Eigen::Index samples_per_word = 32;
}

SampleMask::SampleMask(
    std::span<const Eigen::Index> target_to_source,
    Eigen::Index source_size)
    : source_size_{source_size},
      n_kept_{static_cast<Eigen::Index>(target_to_source.size())}
{
    const auto num_words = static_cast<std::size_t>(
        (source_size + samples_per_word - 1) / samples_per_word);
    words_.assign(num_words, 0);

    for (const Eigen::Index source : target_to_source)
    {
        const auto word = static_cast<std::size_t>(source / samples_per_word);
        const auto slot = static_cast<unsigned>(source % samples_per_word);
        words_[word] |= (std::uint64_t{1} << (2 * slot));
    }
}

}  // namespace gelex::detail
