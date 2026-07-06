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

#ifndef GELEX_BAYES_GENETIC_HALF_NORMAL_PRIOR_STATE_H_
#define GELEX_BAYES_GENETIC_HALF_NORMAL_PRIOR_STATE_H_

#include <array>
#include <string_view>

#include <Eigen/Core>

#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/genetic/prior_state_values.h"
#include "gelex/types/categorical_vector.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::bayes
{

class JointHalfNormalMixtureState final
{
   public:
    static constexpr std::string_view name = "joint_half_normal_mixture";

    JointHalfNormalMixtureState(
        std::array<double, 2> variances,
        const MixtureProportion& proportion,
        const ProbabilityParameter& dominance_positive_probability,
        Eigen::Index num_markers);

    auto variance(GeneticMode mode) -> double&;
    auto variance(GeneticMode mode) const -> const double&;
    auto assignment() -> CategoricalVector& { return mixture_.assignment; }
    auto assignment() const -> const CategoricalVector&
    {
        return mixture_.assignment;
    }
    auto proportion() -> Eigen::VectorXd& { return mixture_.proportion; }
    auto proportion() const -> const Eigen::VectorXd&
    {
        return mixture_.proportion;
    }
    auto dominance_sign() -> DominanceSignState& { return dominance_sign_; }
    auto dominance_sign() const -> const DominanceSignState&
    {
        return dominance_sign_;
    }

    auto visit(FieldVisitor& visitor) -> void;

   private:
    std::array<double, 2> variances_;
    MixtureState mixture_;
    DominanceSignState dominance_sign_;
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_GENETIC_HALF_NORMAL_PRIOR_STATE_H_
