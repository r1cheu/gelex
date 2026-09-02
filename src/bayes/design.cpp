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
#include <array>
#include <cstddef>
#include <fmt/format.h>
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

class GeneticProjection
{
   public:
    GeneticProjection(
        const CompactGenotype& genotype,
        std::span<const gelex::LocusStats> locus_stats,
        const gelex::EncodingSpec& encoding_spec);

    GeneticProjection(const GeneticProjection&) = delete;
    auto operator=(const GeneticProjection&) -> GeneticProjection& = delete;
    GeneticProjection(GeneticProjection&&) noexcept = default;
    auto operator=(GeneticProjection&&) noexcept
        -> GeneticProjection& = default;
    ~GeneticProjection() = default;

    [[nodiscard]] auto rows() const noexcept -> Eigen::Index;
    [[nodiscard]] auto cols() const noexcept -> Eigen::Index;

    [[nodiscard]] auto xtx_diag() const noexcept -> const Eigen::VectorXd&
    {
        return xtx_diag_;
    }

    [[nodiscard]] auto col_var() const noexcept -> const Eigen::RowVectorXd&
    {
        return col_var_;
    }

    [[nodiscard]] auto valid_indices() const noexcept
        -> std::span<const Eigen::Index>
    {
        return valid_indices_;
    }

    [[nodiscard]] auto dot(
        Eigen::Index marker,
        const Eigen::Ref<const Eigen::VectorXd>& values) const noexcept
        -> double;

    auto axpy(
        Eigen::Index marker,
        double scale,
        Eigen::Ref<Eigen::VectorXd> values) const noexcept -> void;

    [[nodiscard]] auto snp_luts() const noexcept -> const gelex::SnpLutMatrix&
    {
        return luts_;
    }
    [[nodiscard]] auto col_covariance(const GeneticProjection& rhs) const
        -> Eigen::RowVectorXd;

   private:
    const CompactGenotype* genotype_;
    gelex::SnpLutMatrix luts_;
    Eigen::VectorXd xtx_diag_;
    Eigen::RowVectorXd col_var_;
    std::vector<Eigen::Index> valid_indices_;
};

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
        const Eigen::Index index = static_cast<Eigen::Index>(marker);
        const gelex::LocusEncoding encoding{
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

class GeneticDesign::Impl
{
   public:
    Impl(gelex::Bed bed, gelex::GenoObserver observer)
        : genotype_(std::move(bed), std::move(observer))
    {
    }

    [[nodiscard]] auto projection() const -> const GeneticProjection&
    {
        if (modes_.size() == 0)
        {
            throw GelexException(
                "GeneticDesign: no genetic projection is available");
        }
        if (modes_.size() > 1)
        {
            throw GelexException(
                "GeneticDesign: multiple genetic projections are available; "
                "specify a mode");
        }
        for (const auto& projection : projections_)
        {
            if (projection)
            {
                return *projection;
            }
        }
        throw GelexException("GeneticDesign: invalid projection state");
    }

    [[nodiscard]] auto projection(GeneticMode mode) const
        -> const GeneticProjection&
    {
        const auto& projection = projections_[mode_index(mode)];
        if (!projection)
        {
            throw GelexException(
                fmt::format(
                    "GeneticDesign: projection for mode {} is not available",
                    mode));
        }
        return *projection;
    }

    CompactGenotype genotype_;
    std::array<std::optional<GeneticProjection>, ALL_GENETIC_MODES.size()>
        projections_;
    GeneticModeSet modes_;
};

GeneticDesign::GeneticDesign(
    gelex::Bed bed,
    GeneticModeSet modes,
    GenotypeMethod geno_method,
    gelex::GenoObserver observer)
    : impl_(std::make_unique<Impl>(std::move(bed), std::move(observer)))
{
    const auto projection_specs
        = encoding_specs_from_method(modes, geno_method);
    for (const auto& spec : projection_specs)
    {
        auto& projection = impl_->projections_[mode_index(spec.effect)];
        projection.emplace(
            impl_->genotype_, impl_->genotype_.locus_stats_, spec);
        impl_->modes_ |= GeneticModeSet{spec.effect};
    }
}

GeneticDesign::GeneticDesign(GeneticDesign&&) noexcept = default;

auto GeneticDesign::operator=(GeneticDesign&&) noexcept
    -> GeneticDesign& = default;

GeneticDesign::~GeneticDesign() = default;

auto GeneticDesign::rows() const noexcept -> Eigen::Index
{
    return impl_->genotype_.rows();
}

auto GeneticDesign::cols() const noexcept -> Eigen::Index
{
    return impl_->genotype_.cols();
}

auto GeneticDesign::modes() const noexcept -> GeneticModeSet
{
    return impl_->modes_;
}

auto GeneticDesign::a1_frequency() const noexcept -> const Eigen::VectorXd&
{
    return impl_->genotype_.a1_frequency();
}

auto GeneticDesign::xtx_diag() const -> const Eigen::VectorXd&
{
    return impl_->projection().xtx_diag();
}

auto GeneticDesign::xtx_diag(GeneticMode mode) const -> const Eigen::VectorXd&
{
    return impl_->projection(mode).xtx_diag();
}

auto GeneticDesign::col_var() const -> const Eigen::RowVectorXd&
{
    return impl_->projection().col_var();
}

auto GeneticDesign::col_var(GeneticMode mode) const -> const Eigen::RowVectorXd&
{
    return impl_->projection(mode).col_var();
}

auto GeneticDesign::valid_indices() const -> std::span<const Eigen::Index>
{
    return impl_->projection().valid_indices();
}

auto GeneticDesign::valid_indices(GeneticMode mode) const
    -> std::span<const Eigen::Index>
{
    return impl_->projection(mode).valid_indices();
}

auto GeneticDesign::dot(
    Eigen::Index marker,
    const Eigen::Ref<const Eigen::VectorXd>& values) const -> double
{
    return impl_->projection().dot(marker, values);
}

auto GeneticDesign::dot(
    GeneticMode mode,
    Eigen::Index marker,
    const Eigen::Ref<const Eigen::VectorXd>& values) const -> double
{
    return impl_->projection(mode).dot(marker, values);
}

auto GeneticDesign::axpy(
    Eigen::Index marker,
    double scale,
    Eigen::Ref<Eigen::VectorXd> values) const -> void
{
    impl_->projection().axpy(marker, scale, values);
}

auto GeneticDesign::axpy(
    GeneticMode mode,
    Eigen::Index marker,
    double scale,
    Eigen::Ref<Eigen::VectorXd> values) const -> void
{
    impl_->projection(mode).axpy(marker, scale, values);
}

auto GeneticDesign::snp_luts() const -> const gelex::SnpLutMatrix&
{
    return impl_->projection().snp_luts();
}

auto GeneticDesign::snp_luts(GeneticMode mode) const
    -> const gelex::SnpLutMatrix&
{
    return impl_->projection(mode).snp_luts();
}

auto GeneticDesign::col_covariance(GeneticMode lhs_mode, GeneticMode rhs_mode)
    const -> Eigen::RowVectorXd
{
    return impl_->projection(lhs_mode).col_covariance(
        impl_->projection(rhs_mode));
}

}  // namespace gelex::bayes
