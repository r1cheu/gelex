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
// for every marker of a BED file, with per-locus code counts derived from them.
class CompactGenotype
{
   public:
    explicit CompactGenotype(
        const gelex::Bed& bed,
        const std::function<void(std::size_t)>& observer = {});

    CompactGenotype(const CompactGenotype&) = delete;
    auto operator=(const CompactGenotype&) -> CompactGenotype& = delete;
    CompactGenotype(CompactGenotype&&) noexcept = default;
    auto operator=(CompactGenotype&&) noexcept -> CompactGenotype& = default;
    ~CompactGenotype() = default;

    [[nodiscard]] auto rows() const noexcept -> Eigen::Index
    {
        return raw_codes_.rows();
    }

    [[nodiscard]] auto cols() const noexcept -> Eigen::Index
    {
        return raw_codes_.cols();
    }

    [[nodiscard]] auto col(Eigen::Index marker) const noexcept
        -> std::span<const std::uint8_t>
    {
        return raw_codes_.col(marker);
    }

    [[nodiscard]] auto locus_stats() const noexcept
        -> std::span<const gelex::LocusStats>
    {
        return locus_stats_;
    }

    [[nodiscard]] auto a1_frequency() const noexcept -> const Eigen::VectorXd&
    {
        return a1_frequency_;
    }

   private:
    using raw_matrix_type = Eigen::
        Array<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

    raw_matrix_type raw_codes_;
    std::vector<gelex::LocusStats> locus_stats_;
    Eigen::VectorXd a1_frequency_;
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_GENOTYPE_COMPACT_GENOTYPE_H_
