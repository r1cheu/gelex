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

#ifndef GELEX_MODEL_BAYES_MODEL_H_
#define GELEX_MODEL_BAYES_MODEL_H_

#include <algorithm>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/states.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

class BayesModel
{
   public:
    BayesModel(
        Eigen::VectorXd phenotype,
        FixedEffect fixed_effects,
        std::vector<bayes::GeneticEffect> genetics,
        bayes::PriorSet priors);

    const FixedEffect& fixed() const { return fixed_; }
    FixedEffect& fixed() { return fixed_; }

    const std::vector<bayes::RandomEffect>& random() const { return random_; }
    std::vector<bayes::RandomEffect>& random() { return random_; }

    const std::vector<bayes::GeneticEffect>& genetics() const
    {
        return genetics_;
    }
    std::vector<bayes::GeneticEffect>& genetics() { return genetics_; }

    const bayes::GeneticEffect* genetic(GeneticKind type) const
    {
        auto it
            = std::ranges::find(genetics_, type, &bayes::GeneticEffect::type);
        return it != genetics_.end() ? &*it : nullptr;
    }
    bayes::GeneticEffect* genetic(GeneticKind type)
    {
        auto it
            = std::ranges::find(genetics_, type, &bayes::GeneticEffect::type);
        return it != genetics_.end() ? &*it : nullptr;
    }

    const bayes::PriorSet& priors() const { return priors_; }
    bayes::PriorSet& priors() { return priors_; }

    const Eigen::VectorXd& phenotype() const { return phenotype_; }

    double phenotype_variance() const { return phenotype_var_; }
    Eigen::Index num_individuals() const { return num_individuals_; }

   private:
    void add_fixed_effect(FixedEffect&& effect);
    void add_random_effect(
        std::string name,
        std::vector<std::string>&& levels,
        Eigen::MatrixXd&& X);

    Eigen::Index num_individuals_{};
    double phenotype_var_{};

    Eigen::VectorXd phenotype_;

    FixedEffect fixed_;
    std::vector<bayes::RandomEffect> random_;
    std::vector<bayes::GeneticEffect> genetics_;
    bayes::PriorSet priors_;
};

class BayesState
{
   public:
    explicit BayesState(const BayesModel& model);

    bayes::FixedState& fixed() { return fixed_; }
    const bayes::FixedState& fixed() const { return fixed_; }
    std::vector<bayes::RandomState>& random() { return random_; }
    const std::vector<bayes::RandomState>& random() const { return random_; }

    const std::vector<bayes::GeneticState>& genetics() const
    {
        return genetics_;
    }
    std::vector<bayes::GeneticState>& genetics() { return genetics_; }

    const bayes::GeneticState* genetic(GeneticKind type) const
    {
        auto it
            = std::ranges::find(genetics_, type, &bayes::GeneticState::type);
        return it != genetics_.end() ? &*it : nullptr;
    }
    bayes::GeneticState* genetic(GeneticKind type)
    {
        auto it
            = std::ranges::find(genetics_, type, &bayes::GeneticState::type);
        return it != genetics_.end() ? &*it : nullptr;
    }

    bayes::ResidualState& residual() { return residual_; }
    const bayes::ResidualState& residual() const { return residual_; }

    void compute_heritability();

   private:
    bayes::FixedState fixed_;
    std::vector<bayes::RandomState> random_;
    std::vector<bayes::GeneticState> genetics_;
    bayes::ResidualState residual_;
};

}  // namespace gelex

#endif  // GELEX_MODEL_BAYES_MODEL_H_
