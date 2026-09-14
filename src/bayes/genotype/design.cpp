// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/bayes/genotype/design.h"

#include <algorithm>
#include <cstddef>
#include <fmt/format.h>
#include <functional>
#include <iterator>
#include <memory>
#include <optional>
#include <string>
#include <utility>

#include "gelex/bayes/genotype/compact_genotype.h"
#include "gelex/bayes/genotype/projection.h"
#include "gelex/data/bed.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/encode/spec.h"
#include "gelex/data/genotype_method.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

namespace gelex::bayes
{

namespace
{

// projections_ is indexed by position in all_genetic_modes, which is the enum
// order; the static_assert pins that assumption.
static_assert(std::ranges::is_sorted(all_genetic_modes));

auto mode_index(GeneticMode mode) -> std::size_t
{
    return static_cast<std::size_t>(std::to_underlying(mode));
}

auto validate_projections(
    const GeneticDesign::projection_array_type& projections,
    const CompactGenotype& genotype) -> void
{
    bool any = false;
    for (const auto& projection : projections)
    {
        if (!projection)
        {
            continue;
        }
        any = true;
        if (&projection->genotype() != &genotype)
        {
            throw GelexException(
                "GeneticDesign: every projection must view the owned genotype");
        }
    }
    if (!any)
    {
        throw GelexException("GeneticDesign: at least one projection required");
    }
}

auto validate_marker_axis(
    const DataFrame<std::string>& marker_metadata,
    const std::optional<MarkerCovariate>& marker_covariate,
    Eigen::Index marker_count) -> void
{
    if (static_cast<Eigen::Index>(marker_metadata.rows()) != marker_count)
    {
        throw GelexException(
            fmt::format(
                "GeneticDesign: marker metadata rows {} != marker count {}",
                marker_metadata.rows(),
                marker_count));
    }
    if (marker_covariate && marker_covariate->X().cols() != marker_count)
    {
        throw GelexException(
            fmt::format(
                "GeneticDesign: marker covariate columns {} != marker count {}",
                marker_covariate->X().cols(),
                marker_count));
    }
}

}  // namespace

GeneticDesign::GeneticDesign(
    std::unique_ptr<CompactGenotype> genotype,
    projection_array_type projections,
    DataFrame<std::string> marker_metadata,
    std::optional<MarkerCovariate> marker_covariate)
    : genotype_{std::move(genotype)},
      marker_metadata_{std::move(marker_metadata)},
      marker_covariate_{std::move(marker_covariate)},
      projections_{std::move(projections)}
{
    if (!genotype_)
    {
        throw GelexException("GeneticDesign: genotype must not be null");
    }
    validate_projections(projections_, *genotype_);
    validate_marker_axis(
        marker_metadata_, marker_covariate_, genotype_->cols());
    if (contains(GeneticMode::A) && contains(GeneticMode::D))
    {
        const auto additive = projection(GeneticMode::A).valid_indices();
        const auto dominance = projection(GeneticMode::D).valid_indices();
        common_valid_indices_.reserve(
            std::min(additive.size(), dominance.size()));
        std::ranges::set_intersection(
            additive, dominance, std::back_inserter(common_valid_indices_));
    }
}

auto GeneticDesign::contains(GeneticMode mode) const -> bool
{
    return projections_.at(mode_index(mode)).has_value();
}

auto GeneticDesign::projection(GeneticMode mode) const
    -> const GeneticProjection&
{
    const auto& value = projections_.at(mode_index(mode));
    if (!value)
    {
        throw GelexException(
            fmt::format(
                "GeneticDesign: projection for mode {} is not available",
                mode));
    }
    return *value;
}

auto GeneticDesign::common_valid_indices() const
    -> std::span<const Eigen::Index>
{
    static_cast<void>(projection(GeneticMode::A));
    static_cast<void>(projection(GeneticMode::D));
    return common_valid_indices_;
}

auto make_genetic_design(
    gelex::Bed bed,
    GeneticModeSet modes,
    GenotypeMethod geno_method,
    std::optional<MarkerCovariate> marker_covariate,
    const std::function<void(std::size_t)>& observer) -> GeneticDesign
{
    auto genotype = std::make_unique<CompactGenotype>(
        make_compact_genotype(bed, observer));
    GeneticDesign::projection_array_type projections;
    for (const GeneticMode mode : modes.each())
    {
        projections.at(mode_index(mode)) = make_genetic_projection(
            *genotype, gelex::encoding_spec_from_method(mode, geno_method));
    }
    return GeneticDesign{
        std::move(genotype),
        std::move(projections),
        std::move(bed).bim(),
        std::move(marker_covariate)};
}

}  // namespace gelex::bayes
