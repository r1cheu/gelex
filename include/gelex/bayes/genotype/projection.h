// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENOTYPE_PROJECTION_H_
#define GELEX_BAYES_GENOTYPE_PROJECTION_H_

#include <Eigen/Core>
#include <span>
#include <vector>

#include "gelex/bayes/genotype/compact_genotype.h"
#include "gelex/data/snp_lut.h"

namespace gelex
{
struct EncodingSpec;
}  // namespace gelex

namespace gelex::bayes
{

// A non-owning encoded view of one CompactGenotype: each raw code column is
// mapped through a per-marker lookup table chosen by the encoding spec. The
// genotype must outlive every projection built on it.
class GeneticProjection
{
   public:
    GeneticProjection(
        const CompactGenotype& genotype,
        const gelex::EncodingSpec& encoding_spec);

    GeneticProjection(const GeneticProjection&) = delete;
    auto operator=(const GeneticProjection&) -> GeneticProjection& = delete;
    GeneticProjection(GeneticProjection&&) noexcept = default;
    auto operator=(GeneticProjection&&) noexcept
        -> GeneticProjection& = default;
    ~GeneticProjection() = default;

    [[nodiscard]] auto rows() const noexcept -> Eigen::Index
    {
        return genotype_->rows();
    }

    [[nodiscard]] auto cols() const noexcept -> Eigen::Index
    {
        return genotype_->cols();
    }

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
        const Eigen::Ref<const Eigen::VectorXd>& rhs) const noexcept -> double;

    auto axpy(
        Eigen::Index marker,
        double scale,
        Eigen::Ref<Eigen::VectorXd> target) const noexcept -> void;

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

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_GENOTYPE_PROJECTION_H_
