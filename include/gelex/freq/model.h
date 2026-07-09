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

#ifndef GELEX_FREQ_MODEL_H_
#define GELEX_FREQ_MODEL_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <span>
#include <vector>

#include "gelex/freq/design.h"
#include "gelex/types/fixed_designs.h"

namespace gelex
{

class FreqModel
{
   public:
    FreqModel(
        Eigen::VectorXd phenotype,
        FixedDesign fixed,
        std::vector<freq::RandomDesign> random);

    auto fixed() const -> const FixedDesign& { return fixed_; }
    auto fixed() -> FixedDesign& { return fixed_; }

    auto random() const -> std::span<const freq::RandomDesign>
    {
        return random_;
    }
    auto random() -> std::span<freq::RandomDesign> { return random_; }

    auto phenotype() const -> const Eigen::VectorXd& { return phenotype_; }
    auto phenotype_variance() const -> double { return phenotype_variance_; }
    auto num_individuals() const -> Eigen::Index { return num_individuals_; }

   private:
    Eigen::Index num_individuals_{};

    Eigen::VectorXd phenotype_;
    double phenotype_variance_;

    FixedDesign fixed_;
    std::vector<freq::RandomDesign> random_;
};

class FreqState
{
   public:
    explicit FreqState(const FreqModel& model);

    freq::FixedState& fixed() { return fixed_; }
    const freq::FixedState& fixed() const { return fixed_; }

    std::span<freq::RandomState> random() { return random_; }
    std::span<const freq::RandomState> random() const { return random_; }

    freq::ResidualState& residual() { return residual_; }
    const freq::ResidualState& residual() const { return residual_; }

    auto Vp() const -> double { return Vp_; }

    void compute_variance_ratio();

   private:
    double Vp_;  // variance of adjusted phenotype
    freq::FixedState fixed_;
    std::vector<freq::RandomState> random_;
    freq::ResidualState residual_;

    void init_variance_components(const FreqModel& model);
};

}  // namespace gelex

#endif  // GELEX_FREQ_MODEL_H_
