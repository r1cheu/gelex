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

#include "gelex/bayes/genotype/projection.h"

#include <Eigen/Core>
#include <cstddef>
#include <ranges>
#include <span>

#include "gelex/data/encode/detail/encoding.h"
#include "gelex/data/encode/stats.h"
#include "gelex/data/encode/types.h"
#include "gelex/exception.h"

#include "bayes/genotype/compact_genotype.h"

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

}  // namespace

namespace gelex::bayes
{

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
        const Eigen::Array4d counts = locus_counts(stats);
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
    return gelex::bayes::dot(genotype_->col(marker), luts_.col(marker), values);
}

auto GeneticProjection::axpy(
    Eigen::Index marker,
    double scale,
    Eigen::Ref<Eigen::VectorXd> values) const noexcept -> void
{
    gelex::bayes::axpy(
        genotype_->col(marker), luts_.col(marker), scale, values);
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
    if (genotype_ != rhs.genotype_)
    {
        throw GelexException(
            "col_covariance: projections must share one compact genotype");
    }

    Eigen::RowVectorXd covariance(cols());
    for (Eigen::Index marker = 0; marker < cols(); ++marker)
    {
        const auto counts = locus_counts(
            genotype_->locus_stats_[static_cast<std::size_t>(marker)]);
        const Eigen::Array4d lhs_values = luts_.col(marker).array();
        const Eigen::Array4d rhs_values = rhs.luts_.col(marker).array();
        const double sample_size = counts.sum();
        const double lhs_mean = (counts * lhs_values).sum() / sample_size;
        const double rhs_mean = (counts * rhs_values).sum() / sample_size;
        covariance[marker]
            = ((counts * lhs_values * rhs_values).sum() / sample_size)
              - (lhs_mean * rhs_mean);
    }
    return covariance;
}

}  // namespace gelex::bayes
