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

#ifndef GELEX_BAYES_GENETIC_PRIOR_STATE_VALUES_H_
#define GELEX_BAYES_GENETIC_PRIOR_STATE_VALUES_H_

#include <string_view>
#include <vector>

#include <Eigen/Core>

#include "gelex/bayes/genetic/parameters.h"

namespace gelex::infra
{
class FieldVisitor;
}

namespace gelex::bayes
{

struct MixtureState
{
    static constexpr std::string_view name = "mixture";

    MixtureState() = default;
    MixtureState(const MixtureProportion& proportion, Eigen::Index num_markers);

    auto visit(infra::FieldVisitor& visitor) -> void;

    Eigen::VectorXi assignment;
    Eigen::VectorXd proportion;
};

struct ComponentState
{
    static constexpr std::string_view name = "component";

    ComponentState() = default;
    ComponentState(Eigen::Index num_components, Eigen::Index num_individuals);

    auto visit(infra::FieldVisitor& visitor) -> void;

    std::vector<Eigen::VectorXd> gebv;
    Eigen::VectorXd gebv_var;
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_GENETIC_PRIOR_STATE_VALUES_H_
