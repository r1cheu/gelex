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

#include "gelex/data/encode/detail/sample_mask.h"

#include <cstddef>

namespace gelex::detail
{

namespace
{
constexpr Eigen::Index SAMPLES_PER_WORD = 32;
}

SampleMask::SampleMask(
    std::span<const Eigen::Index> target_to_source,
    Eigen::Index source_size)
    : source_size_{source_size},
      n_kept_{static_cast<Eigen::Index>(target_to_source.size())}
{
    const auto num_words = static_cast<std::size_t>(
        (source_size + SAMPLES_PER_WORD - 1) / SAMPLES_PER_WORD);
    words_.assign(num_words, 0);

    for (const Eigen::Index source : target_to_source)
    {
        const auto word = static_cast<std::size_t>(source / SAMPLES_PER_WORD);
        const auto slot = static_cast<unsigned>(source % SAMPLES_PER_WORD);
        words_[word] |= (std::uint64_t{1} << (2 * slot));
    }
}

}  // namespace gelex::detail
