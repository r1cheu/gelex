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

#ifndef GELEX_ALGO_INFER_VI_RESULT_H_
#define GELEX_ALGO_INFER_VI_RESULT_H_

#include <algorithm>
#include <vector>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/state.h"
#include "gelex/algo/infer/posterior_summary.h"
#include "gelex/model/bayes/model.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

class BayesModel;

namespace vi
{

struct GeneticSummary
{
    GeneticMode type;
    Eigen::VectorXd coeffs;
};

class Result
{
   public:
    Result(vi::State&& state, const BayesModel& model);

    const FixedSummary& fixed() const { return fixed_; }
    const std::vector<RandomSummary>& random() const { return random_; }

    const std::vector<vi::GeneticSummary>& genetics() const
    {
        return genetics_;
    }
    const vi::GeneticSummary* genetic(GeneticMode type) const
    {
        auto it = std::ranges::find(genetics_, type, &vi::GeneticSummary::type);
        return it != genetics_.end() ? &*it : nullptr;
    }

    const Eigen::VectorXd& allele_freq() const { return p_freq_; }

   private:
    FixedSummary fixed_;
    std::vector<RandomSummary> random_;
    std::vector<vi::GeneticSummary> genetics_;
    Eigen::VectorXd p_freq_;
};

}  // namespace vi

}  // namespace gelex

#endif  // GELEX_ALGO_INFER_VI_RESULT_H_
