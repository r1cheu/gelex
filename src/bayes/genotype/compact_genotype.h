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
