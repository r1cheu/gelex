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

#ifndef GELEX_BAYES_GENOTYPE_DESIGN_H_
#define GELEX_BAYES_GENOTYPE_DESIGN_H_

#include <Eigen/Core>
#include <array>
#include <memory>
#include <optional>
#include <ranges>
#include <span>
#include <vector>

#include "gelex/bayes/genotype/progress.h"
#include "gelex/bayes/genotype/projection.h"
#include "gelex/data/bed.h"
#include "gelex/data/genotype_method.h"
#include "gelex/genetic_mode.h"

namespace gelex::bayes
{

class CompactGenotype;

class GeneticDesign
{
   public:
    GeneticDesign(
        gelex::Bed bed,
        GeneticModeSet modes,
        GenotypeMethod geno_method,
        gelex::GenoObserver observer = {});

    explicit GeneticDesign(gelex::Bed bed, gelex::GenoObserver observer = {});

    GeneticDesign(const GeneticDesign&) = delete;
    auto operator=(const GeneticDesign&) -> GeneticDesign& = delete;
    GeneticDesign(GeneticDesign&&) noexcept;
    auto operator=(GeneticDesign&&) noexcept -> GeneticDesign&;
    ~GeneticDesign();

    [[nodiscard]] auto rows() const noexcept -> Eigen::Index;
    [[nodiscard]] auto cols() const noexcept -> Eigen::Index;

    [[nodiscard]] auto contains(GeneticMode mode) const -> bool;

    [[nodiscard]] auto each_mode() const
    {
        return all_genetic_modes
               | std::views::filter([this](GeneticMode mode)
                                    { return contains(mode); });
    }

    [[nodiscard]] auto a1_frequency() const noexcept -> const Eigen::VectorXd&;

    [[nodiscard]] auto projection(GeneticMode mode) const
        -> const GeneticProjection&;

    [[nodiscard]] auto common_valid_indices() const
        -> std::span<const Eigen::Index>;

   private:
    std::unique_ptr<CompactGenotype> genotype_;
    std::array<std::optional<GeneticProjection>, all_genetic_modes.size()>
        projections_;
    std::vector<Eigen::Index> common_valid_indices_;
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_GENOTYPE_DESIGN_H_
