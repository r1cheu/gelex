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
#include "gelex/model/bayes/vi/states.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

class BayesModel
{
   public:
    BayesModel(
        Eigen::VectorXd phenotype,
        FixedEffect fixed_effects,
        std::vector<bayes::GeneticEffect> genetics);

    const FixedEffect& fixed() const { return fixed_; }

    const std::vector<bayes::RandomEffect>& random() const { return random_; }

    const std::vector<bayes::GeneticEffect>& genetics() const
    {
        return genetics_;
    }

    const bayes::GeneticEffect* genetic(GeneticMode type) const
    {
        auto it
            = std::ranges::find(genetics_, type, &bayes::GeneticEffect::type);
        return it != genetics_.end() ? &*it : nullptr;
    }

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
};

template <typename GeneticStateT>
class InferenceState
{
   public:
    InferenceState(const BayesModel& model, const bayes::Priors& priors);
    InferenceState(
        bayes::FixedState fixed,
        std::vector<bayes::RandomState> random,
        std::vector<GeneticStateT> genetics,
        bayes::ResidualState residual);

    bayes::FixedState& fixed() { return fixed_; }
    const bayes::FixedState& fixed() const { return fixed_; }
    std::vector<bayes::RandomState>& random() { return random_; }
    const std::vector<bayes::RandomState>& random() const { return random_; }

    const std::vector<GeneticStateT>& genetics() const { return genetics_; }
    std::vector<GeneticStateT>& genetics() { return genetics_; }

    const GeneticStateT* genetic(GeneticMode type) const
    {
        auto it = std::ranges::find(genetics_, type, &GeneticStateT::type);
        return it != genetics_.end() ? &*it : nullptr;
    }
    GeneticStateT* genetic(GeneticMode type)
    {
        auto it = std::ranges::find(genetics_, type, &GeneticStateT::type);
        return it != genetics_.end() ? &*it : nullptr;
    }

    bayes::ResidualState& residual() { return residual_; }
    const bayes::ResidualState& residual() const { return residual_; }

    void compute_heritability();

   private:
    bayes::FixedState fixed_;
    std::vector<bayes::RandomState> random_;
    std::vector<GeneticStateT> genetics_;
    bayes::ResidualState residual_;
};

template <typename GeneticStateT>
InferenceState<GeneticStateT>::InferenceState(
    const BayesModel& model,
    const bayes::Priors& priors)
    : fixed_(model.fixed())
{
    for (const auto& effect : model.genetics())
    {
        const auto* prior = priors.genetic(effect.type);
        genetics_.emplace_back(effect, *prior);
    }

    const auto& random_effects = model.random();
    const auto& random_priors = priors.random();
    random_.reserve(random_effects.size());
    for (std::size_t i = 0; i < random_effects.size(); ++i)
    {
        random_.emplace_back(random_effects[i], random_priors[i]);
    }

    residual_.y_adj = model.phenotype().array();
    residual_.variance = priors.residual().init;
}

template <typename GeneticStateT>
InferenceState<GeneticStateT>::InferenceState(
    bayes::FixedState fixed,
    std::vector<bayes::RandomState> random,
    std::vector<GeneticStateT> genetics,
    bayes::ResidualState residual)
    : fixed_(std::move(fixed)),
      random_(std::move(random)),
      genetics_(std::move(genetics)),
      residual_(std::move(residual))
{
}

template <typename GeneticStateT>
void InferenceState<GeneticStateT>::compute_heritability()
{
    double sum_var = 0;

    for (const auto& rand : random_)
    {
        sum_var += rand.variance;
    }

    for (const auto& gen : genetics_)
    {
        sum_var += gen.variance;
    }

    sum_var += residual_.variance;

    for (auto& gen : genetics_)
    {
        gen.heritability = gen.variance / sum_var;
    }
}

namespace mcmc
{
using State = InferenceState<bayes::GeneticState>;
}  // namespace mcmc

namespace vi
{
using State = InferenceState<bayes::vi::GeneticState>;
}  // namespace vi

}  // namespace gelex

#endif  // GELEX_MODEL_BAYES_MODEL_H_
