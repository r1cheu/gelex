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

#include "gelex/bayes/genotype/design.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <fmt/format.h>
#include <iterator>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "gelex/data/bed.h"
#include "gelex/data/encode/spec.h"
#include "gelex/data/genotype_method.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

#include "bayes/genotype/compact_genotype.h"

namespace gelex::bayes
{

namespace
{

[[nodiscard]] auto mode_index(GeneticMode mode) -> std::size_t
{
    switch (mode)
    {
        case GeneticMode::A:
            return 0;
        case GeneticMode::D:
            return 1;
    }
    throw GelexException(fmt::format("Unsupported genetic mode: {}", mode));
}

[[nodiscard]] auto encoding_specs_from_method(
    GeneticModeSet modes,
    GenotypeMethod geno_method) -> std::vector<EncodingSpec>
{
    std::vector<EncodingSpec> specs;
    specs.reserve(modes.size());
    for (const GeneticMode mode : modes.each())
    {
        specs.push_back(gelex::encoding_spec_from_method(mode, geno_method));
    }
    return specs;
}

}  // namespace

GeneticDesign::GeneticDesign(
    gelex::Bed bed,
    GeneticModeSet modes,
    GenotypeMethod geno_method,
    gelex::GenoObserver observer)
    : genotype_{std::make_unique<CompactGenotype>(
          std::move(bed),
          std::move(observer))}
{
    const auto projection_specs
        = encoding_specs_from_method(modes, geno_method);
    for (const auto& spec : projection_specs)
    {
        projections_.at(mode_index(spec.effect))
            = GeneticProjection{*genotype_, genotype_->locus_stats_, spec};
    }
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

GeneticDesign::GeneticDesign(gelex::Bed bed, gelex::GenoObserver observer)
    : genotype_{std::make_unique<CompactGenotype>(
          std::move(bed),
          std::move(observer))}
{
}

GeneticDesign::GeneticDesign(GeneticDesign&&) noexcept = default;

auto GeneticDesign::operator=(GeneticDesign&&) noexcept
    -> GeneticDesign& = default;

GeneticDesign::~GeneticDesign() = default;

auto GeneticDesign::rows() const noexcept -> Eigen::Index
{
    return genotype_->rows();
}

auto GeneticDesign::cols() const noexcept -> Eigen::Index
{
    return genotype_->cols();
}

auto GeneticDesign::contains(GeneticMode mode) const -> bool
{
    return projections_.at(mode_index(mode)).has_value();
}

auto GeneticDesign::a1_frequency() const noexcept -> const Eigen::VectorXd&
{
    return genotype_->a1_frequency();
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

}  // namespace gelex::bayes
