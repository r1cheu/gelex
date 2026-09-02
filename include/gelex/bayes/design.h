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

#ifndef GELEX_BAYES_DESIGN_H_
#define GELEX_BAYES_DESIGN_H_

#include <Eigen/Core>
#include <array>
#include <memory>
#include <optional>
#include <ranges>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "gelex/bayes/genotype_kernel.h"
#include "gelex/data/bed.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/snp_lut.h"
#include "gelex/infra/logging/geno_event.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{
struct EncodingSpec;
class FieldVisitor;
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

struct RandomDesign
{
    RandomDesign(
        std::string name,
        std::vector<std::string> levels,
        Eigen::MatrixXd&& X)
        : name(std::move(name)), levels(std::move(levels)), X(std::move(X))
    {
        XtX_diag = this->X.colwise().squaredNorm();
    }

    std::string name;
    std::vector<std::string> levels;
    Eigen::MatrixXd X;
    Eigen::VectorXd XtX_diag;

    auto visit(FieldVisitor& visitor) const -> void;
};

class GeneticDesign
{
   public:
    GeneticDesign(
        gelex::Bed bed,
        GeneticModeSet modes,
        GenotypeMethod geno_method,
        gelex::GenoObserver observer = {});

    explicit GeneticDesign(gelex::Bed bed, gelex::GenoObserver observer = {});

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

    [[nodiscard]] auto projection(GeneticMode mode) const
        -> const GeneticProjection&;

    [[nodiscard]] auto common_valid_indices() const
        -> std::span<const Eigen::Index>;

   private:
    std::unique_ptr<CompactGenotype> genotype_;
    std::array<std::optional<GeneticProjection>, all_genetic_modes.size()>
        projections_;
    std::vector<Eigen::Index> common_valid_indices_;
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_DESIGN_H_
