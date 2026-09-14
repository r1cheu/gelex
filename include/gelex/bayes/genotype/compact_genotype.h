// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENOTYPE_COMPACT_GENOTYPE_H_
#define GELEX_BAYES_GENOTYPE_COMPACT_GENOTYPE_H_

#include <Eigen/Core>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <span>
#include <vector>

#include "gelex/data/bed.h"
#include "gelex/data/encode/stats.h"

namespace gelex::bayes
{

// Column-major raw genotype codes (0 = A1A1, 1 = missing, 2 = A1A2, 3 = A2A2)
// for every marker, with per-locus code counts and A1 frequencies indexed by
// the same marker axis.
class CompactGenotype
{
   public:
    using raw_matrix_type = Eigen::
        Array<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

    // Throws GelexException when the three inputs disagree on the marker count.
    CompactGenotype(
        raw_matrix_type raw_codes,
        std::vector<gelex::LocusStats> locus_stats,
        Eigen::VectorXd a1_frequency);

    CompactGenotype(const CompactGenotype&) = delete;
    auto operator=(const CompactGenotype&) -> CompactGenotype& = delete;
    CompactGenotype(CompactGenotype&&) noexcept = default;
    auto operator=(CompactGenotype&&) noexcept -> CompactGenotype& = default;
    ~CompactGenotype() = default;

    auto rows() const noexcept -> Eigen::Index { return raw_codes_.rows(); }

    auto cols() const noexcept -> Eigen::Index { return raw_codes_.cols(); }

    auto col(Eigen::Index marker) const noexcept
        -> std::span<const std::uint8_t>
    {
        return raw_codes_.col(marker);
    }

    auto locus_stats() const noexcept -> std::span<const gelex::LocusStats>
    {
        return locus_stats_;
    }

    auto a1_frequency() const noexcept -> const Eigen::VectorXd&
    {
        return a1_frequency_;
    }

   private:
    raw_matrix_type raw_codes_;
    std::vector<gelex::LocusStats> locus_stats_;
    Eigen::VectorXd a1_frequency_;
};

// Decodes every marker of the BED's current sample selection; observer is
// notified with the number of markers completed so far.
[[nodiscard]] auto make_compact_genotype(
    const gelex::Bed& bed,
    const std::function<void(std::size_t)>& observer = {}) -> CompactGenotype;

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_GENOTYPE_COMPACT_GENOTYPE_H_
