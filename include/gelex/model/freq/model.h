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

#ifndef GELEX_MODEL_FREQ_FREQMODEL_H_
#define GELEX_MODEL_FREQ_FREQMODEL_H_

#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "../src/types/fixed_effects.h"
#include "../src/types/freq_effect.h"
#include "gelex/data/data_pipe.h"

namespace gelex
{

class FreqModel
{
   public:
    explicit FreqModel(DataPipe& data_pipe);
    auto fixed() const -> const FixedEffect& { return fixed_; }
    auto fixed() -> FixedEffect& { return fixed_; }

    auto random() const -> const std::vector<freq::RandomEffect>&
    {
        return random_;
    }
    auto random() -> std::vector<freq::RandomEffect>& { return random_; }

    auto genetic() const -> const std::vector<freq::GeneticEffect>&
    {
        return genetic_;
    }
    auto genetic() -> std::vector<freq::GeneticEffect>& { return genetic_; }

    auto phenotype() const -> const Eigen::VectorXd& { return phenotype_; }
    auto phenotype_variance() const -> double { return phenotype_variance_; }
    auto num_individuals() const -> Eigen::Index { return num_individuals_; }

   private:
    Eigen::Index num_individuals_{};

    Eigen::VectorXd phenotype_;
    double phenotype_variance_;

    FixedEffect fixed_;
    std::vector<freq::RandomEffect> random_;
    std::vector<freq::GeneticEffect> genetic_;
};

class FreqState
{
   public:
    explicit FreqState(const FreqModel& model);

    freq::FixedState& fixed() { return fixed_; }
    const freq::FixedState& fixed() const { return fixed_; }

    std::vector<freq::RandomState>& random() { return random_; }
    const std::vector<freq::RandomState>& random() const { return random_; }

    std::vector<freq::GeneticState>& genetic() { return genetic_; }
    const std::vector<freq::GeneticState>& genetic() const { return genetic_; }

    freq::ResidualState& residual() { return residual_; }
    const freq::ResidualState& residual() const { return residual_; }

    void compute_heritability();

   private:
    double phenotype_variance_;
    freq::FixedState fixed_;
    std::vector<freq::RandomState> random_;
    std::vector<freq::GeneticState> genetic_;
    freq::ResidualState residual_;

    void init_variance_components(const FreqModel& model);
};

}  // namespace gelex

#endif  // GELEX_MODEL_FREQ_FREQMODEL_H_
