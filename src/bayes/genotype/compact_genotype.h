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

class GeneticDesign;
class GeneticProjection;

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

    [[nodiscard]] auto rows() const noexcept -> Eigen::Index;
    [[nodiscard]] auto cols() const noexcept -> Eigen::Index;
    [[nodiscard]] auto size_bytes() const noexcept -> std::size_t;

    [[nodiscard]] auto a1_frequency() const noexcept -> const Eigen::VectorXd&
    {
        return a1_frequency_;
    }

   private:
    using RawMatrix = Eigen::
        Array<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

    [[nodiscard]] auto col(Eigen::Index index) const noexcept
        -> std::span<const std::uint8_t>;

    RawMatrix raw_codes_;
    std::vector<gelex::LocusStats> locus_stats_;
    Eigen::VectorXd a1_frequency_;

    friend class GeneticDesign;
    friend class GeneticProjection;
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_GENOTYPE_COMPACT_GENOTYPE_H_
