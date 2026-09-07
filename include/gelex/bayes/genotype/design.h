// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENOTYPE_DESIGN_H_
#define GELEX_BAYES_GENOTYPE_DESIGN_H_

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <ranges>
#include <span>
#include <string>
#include <vector>

#include "gelex/bayes/genetic/marker_covariate.h"
#include "gelex/bayes/genotype/projection.h"
#include "gelex/data/bed.h"
#include "gelex/data/dataframe/dataframe.h"
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
        std::optional<MarkerCovariate> marker_covariate = std::nullopt,
        const std::function<void(std::size_t)>& observer = {});

    explicit GeneticDesign(
        gelex::Bed bed,
        std::optional<MarkerCovariate> marker_covariate = std::nullopt,
        const std::function<void(std::size_t)>& observer = {});

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

    [[nodiscard]] auto marker_metadata() const noexcept
        -> const DataFrame<std::string>&
    {
        return marker_metadata_;
    }

    [[nodiscard]] auto marker_covariate() const noexcept
        -> const std::optional<MarkerCovariate>&
    {
        return marker_covariate_;
    }

    [[nodiscard]] auto projection(GeneticMode mode) const
        -> const GeneticProjection&;

    [[nodiscard]] auto common_valid_indices() const
        -> std::span<const Eigen::Index>;

   private:
    std::unique_ptr<CompactGenotype> genotype_;
    DataFrame<std::string> marker_metadata_;
    std::optional<MarkerCovariate> marker_covariate_;
    std::array<std::optional<GeneticProjection>, all_genetic_modes.size()>
        projections_;
    std::vector<Eigen::Index> common_valid_indices_;
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_GENOTYPE_DESIGN_H_
