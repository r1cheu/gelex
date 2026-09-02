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
#include <memory>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "gelex/data/bed.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/snp_lut.h"
#include "gelex/infra/logging/geno_event.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{
class FieldVisitor;
}  // namespace gelex

namespace gelex::bayes
{

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

    GeneticDesign(const GeneticDesign&) = delete;
    auto operator=(const GeneticDesign&) -> GeneticDesign& = delete;
    GeneticDesign(GeneticDesign&&) noexcept;
    auto operator=(GeneticDesign&&) noexcept -> GeneticDesign&;
    ~GeneticDesign();

    [[nodiscard]] auto rows() const noexcept -> Eigen::Index;
    [[nodiscard]] auto cols() const noexcept -> Eigen::Index;
    [[nodiscard]] auto modes() const noexcept -> GeneticModeSet;

    [[nodiscard]] auto a1_frequency() const noexcept -> const Eigen::VectorXd&;

    [[nodiscard]] auto xtx_diag() const -> const Eigen::VectorXd&;
    [[nodiscard]] auto xtx_diag(GeneticMode mode) const
        -> const Eigen::VectorXd&;

    [[nodiscard]] auto col_var() const -> const Eigen::RowVectorXd&;
    [[nodiscard]] auto col_var(GeneticMode mode) const
        -> const Eigen::RowVectorXd&;

    [[nodiscard]] auto valid_indices() const -> std::span<const Eigen::Index>;
    [[nodiscard]] auto valid_indices(GeneticMode mode) const
        -> std::span<const Eigen::Index>;

    [[nodiscard]] auto dot(
        Eigen::Index marker,
        const Eigen::Ref<const Eigen::VectorXd>& values) const -> double;
    [[nodiscard]] auto dot(
        GeneticMode mode,
        Eigen::Index marker,
        const Eigen::Ref<const Eigen::VectorXd>& values) const -> double;

    auto axpy(
        Eigen::Index marker,
        double scale,
        Eigen::Ref<Eigen::VectorXd> values) const -> void;
    auto axpy(
        GeneticMode mode,
        Eigen::Index marker,
        double scale,
        Eigen::Ref<Eigen::VectorXd> values) const -> void;

    [[nodiscard]] auto snp_luts() const -> const gelex::SnpLutMatrix&;
    [[nodiscard]] auto snp_luts(GeneticMode mode) const
        -> const gelex::SnpLutMatrix&;

    [[nodiscard]] auto col_covariance(
        GeneticMode lhs_mode,
        GeneticMode rhs_mode) const -> Eigen::RowVectorXd;

   private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_DESIGN_H_
