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

#ifndef GELEX_BAYES_GENETIC_HALF_NORMAL_PRIOR_H_
#define GELEX_BAYES_GENETIC_HALF_NORMAL_PRIOR_H_

#include <string_view>

#include <Eigen/Core>

#include "gelex/bayes/genetic/half_normal_prior_state.h"
#include "gelex/bayes/genetic/parameters.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::bayes
{

class JointHalfNormalMixturePrior final
{
   public:
    static constexpr std::string_view name = "joint_half_normal_mixture";

    JointHalfNormalMixturePrior(
        JointSharedMarkerVariance variance,
        MixtureProportion proportion,
        ProbabilityParameter dominance_positive_probability);

    auto variance(GeneticMode mode) -> SharedMarkerVariance&;
    auto variance(GeneticMode mode) const -> const SharedMarkerVariance&;
    auto variances() -> JointSharedMarkerVariance& { return variances_; }
    auto variances() const -> const JointSharedMarkerVariance&
    {
        return variances_;
    }
    auto proportion() -> MixtureProportion& { return proportion_; }
    auto proportion() const -> const MixtureProportion& { return proportion_; }
    auto dominance_positive_probability() -> ProbabilityParameter&
    {
        return dominance_positive_probability_;
    }
    auto dominance_positive_probability() const -> const ProbabilityParameter&
    {
        return dominance_positive_probability_;
    }

    auto visit(FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> JointHalfNormalMixtureState;

   private:
    JointSharedMarkerVariance variances_;
    MixtureProportion proportion_;
    ProbabilityParameter dominance_positive_probability_;
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_GENETIC_HALF_NORMAL_PRIOR_H_
