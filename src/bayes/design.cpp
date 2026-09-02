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

#include "gelex/bayes/design.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cstddef>
#include <fmt/format.h>
#include <iterator>
#include <memory>
#include <optional>
#include <ranges>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "gelex/bayes/genotype_kernel.h"
#include "gelex/data/bed.h"
#include "gelex/data/encode/detail/encoding.h"
#include "gelex/data/encode/spec.h"
#include "gelex/data/encode/stats.h"
#include "gelex/data/encode/types.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/snp_lut.h"
#include "gelex/exception.h"
#include "gelex/infra/field_flag.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/infra/logging/geno_event.h"
#include "gelex/types/genetic_mode.h"

#include "bayes/compact_genotype.h"

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

auto RandomDesign::visit(FieldVisitor& visitor) const -> void
{
    std::vector<std::string> coefficient_names;
    coefficient_names.reserve(levels.size());
    for (const auto& level : levels)
    {
        coefficient_names.push_back(name.empty() ? level : name + "_" + level);
    }

    visitor.on(
        "coeffs_names",
        std::span<const std::string>{coefficient_names},
        FieldFlag::report);
    visitor.on("variance_name", fmt::format("σ²_{}", name), FieldFlag::report);
}

GeneticProjection::GeneticProjection(
    const CompactGenotype& genotype,
    std::span<const gelex::LocusStats> locus_stats,
    const gelex::EncodingSpec& encoding_spec)
    : genotype_(&genotype),
      luts_(4, genotype.cols()),
      xtx_diag_(genotype.cols()),
      col_var_(genotype.cols())
{
    luts_.setZero();
    xtx_diag_.setZero();
    col_var_.setZero();
    valid_indices_.reserve(static_cast<std::size_t>(genotype.cols()));

    for (const auto [marker, stats] : std::views::enumerate(locus_stats))
    {
        const auto index = static_cast<Eigen::Index>(marker);
        const auto encoding = gelex::LocusEncoding{
            gelex::detail::make_locus_encoding(index, stats, encoding_spec)};
        if (!encoding.valid)
        {
            continue;
        }

        luts_.col(index) = encoding.lut;
        const Eigen::Array4d counts{
            static_cast<double>(stats.nA1A1),
            static_cast<double>(stats.n_missing),
            static_cast<double>(stats.nA1A2),
            static_cast<double>(stats.nA2A2)};
        const auto values = luts_.col(index);
        const double sum = (counts * values).sum();
        const double sum_sq = (counts * values.square()).sum();
        const double sample_size = counts.sum();
        xtx_diag_[index] = sum_sq;
        col_var_[index] = (sum_sq / sample_size)
                          - ((sum / sample_size) * (sum / sample_size));
        valid_indices_.push_back(index);
    }
}

auto GeneticProjection::rows() const noexcept -> Eigen::Index
{
    return genotype_->rows();
}

auto GeneticProjection::cols() const noexcept -> Eigen::Index
{
    return genotype_->cols();
}

auto GeneticProjection::dot(
    Eigen::Index marker,
    const Eigen::Ref<const Eigen::VectorXd>& values) const noexcept -> double
{
    return gelex::bayes::dot(
        genotype_->col(marker),
        luts_.col(marker),
        std::span<const double>{
            values.data(), static_cast<std::size_t>(values.size())});
}

auto GeneticProjection::axpy(
    Eigen::Index marker,
    double scale,
    Eigen::Ref<Eigen::VectorXd> values) const noexcept -> void
{
    gelex::bayes::axpy(
        genotype_->col(marker),
        luts_.col(marker),
        scale,
        std::span<double>{
            values.data(), static_cast<std::size_t>(values.size())});
}

auto GeneticProjection::axpy(
    Eigen::Index marker,
    std::span<const AxpyTarget> targets) const noexcept -> void
{
    gelex::bayes::axpy(genotype_->col(marker), luts_.col(marker), targets);
}

auto GeneticProjection::col_covariance(const GeneticProjection& rhs) const
    -> Eigen::RowVectorXd
{
    if (rows() == 0)
    {
        throw GelexException("col_covariance: projections have no rows");
    }

    Eigen::RowVectorXd covariance(cols());
    Eigen::VectorXd lhs_column(rows());
    Eigen::VectorXd rhs_column(rhs.rows());
    for (Eigen::Index marker = 0; marker < cols(); ++marker)
    {
        lhs_column.setZero();
        rhs_column.setZero();
        axpy(marker, 1.0, lhs_column);
        rhs.axpy(marker, 1.0, rhs_column);
        covariance[marker] = (lhs_column.array() * rhs_column.array()).mean()
                             - (lhs_column.mean() * rhs_column.mean());
    }
    return covariance;
}

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
