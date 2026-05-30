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

#ifndef GELEX_MODEL_BAYES_PRIOR_STATE_VALUES_H_
#define GELEX_MODEL_BAYES_PRIOR_STATE_VALUES_H_

#include <string_view>
#include <vector>

#include <Eigen/Core>

#include "gelex/model/bayes/genetic_prior_parameters.h"

namespace gelex::infra
{
class FieldVisitor;
}

namespace gelex::bayes
{

struct AssignmentState
{
    static constexpr std::string_view name = "mixture_assignment";

    AssignmentState() = default;
    AssignmentState(Eigen::Index num_components, Eigen::Index num_markers);

    auto visit(infra::FieldVisitor& visitor) -> void;

    Eigen::VectorXi assignment;
    Eigen::VectorXi count;
};

struct SampledProportionState
{
    static constexpr std::string_view name = "sampled_mixture_proportion";

    SampledProportionState() = default;
    SampledProportionState(
        const SampledProportion& proportion,
        Eigen::Index num_markers);

    auto visit(infra::FieldVisitor& visitor) -> void;

    AssignmentState assignment;
    Eigen::VectorXd value;
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

#endif  // GELEX_MODEL_BAYES_PRIOR_STATE_VALUES_H_
