// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/bayes/genotype/projection.h"

#include <Eigen/Core>
#include <algorithm>
#include <cstddef>
#include <fmt/format.h>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

#include "gelex/bayes/genotype/compact_genotype.h"
#include "gelex/bayes/genotype/operations.h"
#include "gelex/data/encode/detail/encoding.h"
#include "gelex/data/encode/stats.h"
#include "gelex/data/encode/types.h"
#include "gelex/data/snp_lut.h"
#include "gelex/exception.h"

namespace
{

auto locus_counts(const gelex::LocusStats& stats) -> Eigen::Array4d
{
    return {
        static_cast<double>(stats.nA1A1),
        static_cast<double>(stats.n_missing),
        static_cast<double>(stats.nA1A2),
        static_cast<double>(stats.nA2A2)};
}

auto validate_valid_indices(
    std::span<const Eigen::Index> valid_indices,
    Eigen::Index marker_count) -> void
{
    const bool in_range = std::ranges::all_of(
        valid_indices,
        [marker_count](Eigen::Index index)
        { return index >= 0 && index < marker_count; });
    if (!in_range)
    {
        throw gelex::GelexException(
            fmt::format(
                "GeneticProjection: valid index out of range for {} markers",
                marker_count));
    }
    if (std::ranges::adjacent_find(valid_indices, std::ranges::greater_equal{})
        != valid_indices.end())
    {
        throw gelex::GelexException(
            "GeneticProjection: valid indices must be strictly increasing");
    }
}

}  // namespace

namespace gelex::bayes
{

GeneticProjection::GeneticProjection(
    const CompactGenotype& genotype,
    gelex::SnpLutMatrix luts,
    std::vector<Eigen::Index> valid_indices)
    : genotype_(&genotype),
      luts_(std::move(luts)),
      xtx_diag_(Eigen::VectorXd::Zero(genotype.cols())),
      col_var_(Eigen::RowVectorXd::Zero(genotype.cols())),
      valid_indices_(std::move(valid_indices))
{
    if (luts_.cols() != genotype.cols())
    {
        throw GelexException(
            fmt::format(
                "GeneticProjection: lookup table has {} columns but genotype "
                "has {} markers",
                luts_.cols(),
                genotype.cols()));
    }
    validate_valid_indices(valid_indices_, genotype.cols());

    const auto locus_stats = genotype.locus_stats();
    for (const Eigen::Index marker : valid_indices_)
    {
        const Eigen::Array4d counts
            = locus_counts(locus_stats[static_cast<std::size_t>(marker)]);
        const auto values = luts_.col(marker);
        const double sum = (counts * values).sum();
        const double sum_sq = (counts * values.square()).sum();
        const double sample_size = counts.sum();
        xtx_diag_[marker] = sum_sq;
        col_var_[marker] = (sum_sq / sample_size)
                           - ((sum / sample_size) * (sum / sample_size));
    }
}

auto GeneticProjection::dot(
    Eigen::Index marker,
    const Eigen::Ref<const Eigen::VectorXd>& rhs) const noexcept -> double
{
    return gelex::bayes::dot(genotype_->col(marker), luts_.col(marker), rhs);
}

auto GeneticProjection::axpy(
    Eigen::Index marker,
    double scale,
    Eigen::Ref<Eigen::VectorXd> target) const noexcept -> void
{
    gelex::bayes::axpy(
        genotype_->col(marker), luts_.col(marker), scale, target);
}

auto GeneticProjection::col_covariance(const GeneticProjection& rhs) const
    -> Eigen::RowVectorXd
{
    if (genotype_ != rhs.genotype_)
    {
        throw GelexException(
            "col_covariance: projections must share one compact genotype");
    }

    const auto locus_stats = genotype_->locus_stats();
    Eigen::RowVectorXd covariance(cols());
    for (const auto [marker, stats] : std::views::enumerate(locus_stats))
    {
        const auto index = static_cast<Eigen::Index>(marker);
        const auto counts = locus_counts(stats);
        const Eigen::Array4d lhs_values = luts_.col(index).array();
        const Eigen::Array4d rhs_values = rhs.luts_.col(index).array();
        const double sample_size = counts.sum();
        const double lhs_mean = (counts * lhs_values).sum() / sample_size;
        const double rhs_mean = (counts * rhs_values).sum() / sample_size;
        covariance[index]
            = ((counts * lhs_values * rhs_values).sum() / sample_size)
              - (lhs_mean * rhs_mean);
    }
    return covariance;
}

auto make_genetic_projection(
    const CompactGenotype& genotype,
    const gelex::EncodingSpec& spec) -> GeneticProjection
{
    gelex::SnpLutMatrix luts = gelex::SnpLutMatrix::Zero(4, genotype.cols());
    std::vector<Eigen::Index> valid_indices;
    valid_indices.reserve(static_cast<std::size_t>(genotype.cols()));

    for (const auto [marker, stats] :
         std::views::enumerate(genotype.locus_stats()))
    {
        const auto index = static_cast<Eigen::Index>(marker);
        const auto encoding = gelex::LocusEncoding{
            gelex::detail::make_locus_encoding(index, stats, spec)};
        if (!encoding.valid)
        {
            continue;
        }
        luts.col(index) = encoding.lut;
        valid_indices.push_back(index);
    }
    return GeneticProjection{
        genotype, std::move(luts), std::move(valid_indices)};
}

}  // namespace gelex::bayes
