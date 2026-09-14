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
#include "gelex/bayes/genotype/compact_genotype.h"
#include "gelex/bayes/genotype/projection.h"
#include "gelex/data/bed.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/snp_lut.h"
#include "gelex/genetic_mode.h"

namespace gelex::bayes
{

// Owns one CompactGenotype together with the projection of every requested
// genetic mode onto it, plus the marker-aligned metadata and covariates.
class GeneticDesign
{
   public:
    // Indexed by the underlying value of GeneticMode; empty slots are modes
    // the design does not model.
    using projection_array_type = std::
        array<std::optional<GeneticProjection>, all_genetic_modes.size()>;

    // Every projection must view *genotype, at least one must be present, and
    // marker_metadata / marker_covariate must span genotype->cols() markers.
    // Throws GelexException otherwise.
    GeneticDesign(
        std::unique_ptr<CompactGenotype> genotype,
        projection_array_type projections,
        DataFrame<std::string> marker_metadata,
        std::optional<MarkerCovariate> marker_covariate = std::nullopt);

    GeneticDesign(const GeneticDesign&) = delete;
    auto operator=(const GeneticDesign&) -> GeneticDesign& = delete;
    GeneticDesign(GeneticDesign&&) noexcept = default;
    auto operator=(GeneticDesign&&) noexcept -> GeneticDesign& = default;
    ~GeneticDesign() = default;

    auto rows() const noexcept -> Eigen::Index { return genotype_->rows(); }

    auto cols() const noexcept -> Eigen::Index { return genotype_->cols(); }

    [[nodiscard]] auto contains(GeneticMode mode) const -> bool;

    auto each_mode() const
    {
        return all_genetic_modes
               | std::views::filter([this](GeneticMode mode)
                                    { return contains(mode); });
    }

    auto a1_frequency() const noexcept -> const Eigen::VectorXd&
    {
        return genotype_->a1_frequency();
    }

    auto marker_metadata() const noexcept -> const DataFrame<std::string>&
    {
        return marker_metadata_;
    }

    auto marker_covariate() const noexcept
        -> const std::optional<MarkerCovariate>&
    {
        return marker_covariate_;
    }

    [[nodiscard]] auto projection(GeneticMode mode) const
        -> const GeneticProjection&;

    [[nodiscard]] auto common_valid_indices() const
        -> std::span<const Eigen::Index>;

   private:
    // Heap-allocated so its address survives moves of the design; every
    // projection below holds a non-owning pointer to it.
    std::unique_ptr<CompactGenotype> genotype_;
    DataFrame<std::string> marker_metadata_;
    std::optional<MarkerCovariate> marker_covariate_;
    projection_array_type projections_;
    std::vector<Eigen::Index> common_valid_indices_;
};

// Decodes the BED's current sample selection, takes its bim as marker
// metadata, and projects every mode in modes under geno_method. observer is
// notified with the number of markers decoded so far. The rvalue overload
// moves the bim out of the BED; the lvalue overload clones it.
[[nodiscard]] auto make_genetic_design(
    gelex::Bed&& bed,
    GeneticModeSet modes,
    GenotypeMethod geno_method,
    std::optional<MarkerCovariate> marker_covariate = std::nullopt,
    const std::function<void(std::size_t)>& observer = {}) -> GeneticDesign;

[[nodiscard]] auto make_genetic_design(
    const gelex::Bed& bed,
    GeneticModeSet modes,
    GenotypeMethod geno_method,
    std::optional<MarkerCovariate> marker_covariate = std::nullopt,
    const std::function<void(std::size_t)>& observer = {}) -> GeneticDesign;

// Rebuilds the design of a trained model from its saved lookup tables (one
// per mode, 4 x markers). The BED must carry the training markers in the same
// order and allele orientation; only the marker count can be verified here.
[[nodiscard]] auto make_genetic_design(
    const gelex::Bed& bed,
    const ModeMap<gelex::SnpLutMatrix>& luts,
    std::optional<MarkerCovariate> marker_covariate = std::nullopt,
    const std::function<void(std::size_t)>& observer = {}) -> GeneticDesign;

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_GENOTYPE_DESIGN_H_
