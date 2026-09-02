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

#include "bayes/compact_genotype.h"

#include <Eigen/Core>
#include <cstddef>
#include <cstdint>
#include <span>

#include "gelex/data/encode/encoder.h"
#include "gelex/infra/logging/notify.h"

namespace gelex::bayes
{

CompactGenotype::CompactGenotype(gelex::Bed bed, gelex::GenoObserver observer)
    : raw_codes_(bed.num_samples(), bed.num_snps()),
      locus_stats_(static_cast<std::size_t>(bed.num_snps())),
      a1_frequency_(bed.num_snps())
{
    const gelex::LocusEncoder encoder{bed};
    for (Eigen::Index marker = 0; marker < cols(); ++marker)
    {
        auto* column = raw_codes_.col(marker).data();
        encoder.decode_raw_codes(
            marker,
            std::span<std::uint8_t>{column, static_cast<std::size_t>(rows())});

        auto& stats = locus_stats_[static_cast<std::size_t>(marker)];
        for (Eigen::Index sample = 0; sample < rows(); ++sample)
        {
            switch (raw_codes_(sample, marker))
            {
                case std::uint8_t{0}:
                    ++stats.nA1A1;
                    break;
                case std::uint8_t{1}:
                    ++stats.n_missing;
                    break;
                case std::uint8_t{2}:
                    ++stats.nA1A2;
                    break;
                case std::uint8_t{3}:
                    ++stats.nA2A2;
                    break;
            }
        }
        a1_frequency_[marker] = stats.has_nonmissing() ? stats.A1freq() : 0.0;
        notify(
            observer,
            gelex::GenotypeProgressEvent{
                static_cast<std::size_t>(marker + 1),
                static_cast<std::size_t>(cols()),
                false});
    }
    notify(
        observer,
        gelex::GenotypeProgressEvent{
            static_cast<std::size_t>(cols()),
            static_cast<std::size_t>(cols()),
            true});
}

auto CompactGenotype::rows() const noexcept -> Eigen::Index
{
    return raw_codes_.rows();
}

auto CompactGenotype::cols() const noexcept -> Eigen::Index
{
    return raw_codes_.cols();
}

auto CompactGenotype::size_bytes() const noexcept -> std::size_t
{
    return static_cast<std::size_t>(raw_codes_.size());
}

auto CompactGenotype::col(Eigen::Index index) const noexcept
    -> std::span<const std::uint8_t>
{
    return std::span<const std::uint8_t>{
        raw_codes_.col(index).data(), static_cast<std::size_t>(rows())};
}

}  // namespace gelex::bayes
