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

#ifndef GELEX_BAYES_GENETIC_PROJECTION_H_
#define GELEX_BAYES_GENETIC_PROJECTION_H_

#include <Eigen/Core>
#include <span>
#include <vector>

#include "gelex/bayes/genotype_operations.h"
#include "gelex/data/snp_lut.h"

namespace gelex
{
struct EncodingSpec;
struct LocusStats;
}  // namespace gelex

namespace gelex::bayes
{

class CompactGenotype;

class GeneticProjection
{
   public:
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

    auto axpy(Eigen::Index marker, std::span<const AxpyTarget> targets)
        const noexcept -> void;

    [[nodiscard]] auto snp_luts() const noexcept -> const gelex::SnpLutMatrix&
    {
        return luts_;
    }

    [[nodiscard]] auto col_covariance(const GeneticProjection& rhs) const
        -> Eigen::RowVectorXd;

   private:
    GeneticProjection(
        const CompactGenotype& genotype,
        std::span<const gelex::LocusStats> locus_stats,
        const gelex::EncodingSpec& encoding_spec);

    const CompactGenotype* genotype_;
    gelex::SnpLutMatrix luts_;
    Eigen::VectorXd xtx_diag_;
    Eigen::RowVectorXd col_var_;
    std::vector<Eigen::Index> valid_indices_;

    friend class GeneticDesign;
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_GENETIC_PROJECTION_H_
