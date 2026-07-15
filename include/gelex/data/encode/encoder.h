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

#ifndef GELEX_DATA_ENCODE_ENCODER_H_
#define GELEX_DATA_ENCODE_ENCODER_H_

#include <Eigen/Core>
#include <cassert>
#include <cstddef>
#include <span>

#include "gelex/data/bed.h"
#include "gelex/data/detail/bed_source.h"
#include "gelex/data/detail/index_projection.h"
#include "gelex/data/encode/detail/count.h"
#include "gelex/data/encode/detail/encoding.h"
#include "gelex/data/encode/detail/expand.h"
#include "gelex/data/encode/detail/sample_mask.h"
#include "gelex/data/encode/stats.h"
#include "gelex/data/encode/types.h"

namespace gelex
{

// Counts, encodes, and expands any variant straight from a Bed's 2-bit-packed
// source, never materializing a dosage matrix. Every method is const and
// side-effect-free per variant, so distinct variants may be processed
// concurrently. count() is split from encoding() so one count can feed several
// specs' encodings.
class LocusEncoder
{
   public:
    explicit LocusEncoder(const Bed& bed)
        : source_{bed.bed_source_},
          mask_{bed.index_projection_},
          target_to_source_{bed.index_projection_.target_to_source()}
    {
    }

    [[nodiscard]] auto target_size() const noexcept -> Eigen::Index
    {
        return static_cast<Eigen::Index>(target_to_source_.size());
    }

    [[nodiscard]] auto count(Eigen::Index variant) const -> LocusStats
    {
        return detail::count_genotypes(
            source_[static_cast<std::size_t>(variant)], mask_);
    }

    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    [[nodiscard]] auto encoding(
        Eigen::Index variant,
        const LocusStats& stats,
        const EncodingSpec& spec) const -> LocusEncoding
    {
        return detail::make_locus_encoding(variant, stats, spec);
    }

    // encoding may come from encoding() on this variant, or externally (e.g.
    // training codes for prediction, with any allele flip baked into
    // encoding.code); stats are not recomputed from the current samples.
    auto expand(
        Eigen::Index variant,
        const LocusEncoding& encoding,
        Eigen::Ref<Eigen::VectorXd> out) const -> void
    {
        assert(out.size() == target_size());
        detail::expand_encoded_column(
            source_[static_cast<std::size_t>(variant)],
            target_to_source_,
            encoding,
            out);
    }

   private:
    const detail::BedSource& source_;
    detail::SampleMask mask_;
    std::span<const Eigen::Index> target_to_source_;
};

}  // namespace gelex

#endif  // GELEX_DATA_ENCODE_ENCODER_H_
