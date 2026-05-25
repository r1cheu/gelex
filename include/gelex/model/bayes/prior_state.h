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

#ifndef GELEX_MODEL_BAYES_PRIOR_STATE_H_
#define GELEX_MODEL_BAYES_PRIOR_STATE_H_

#include <cstddef>

#include <Eigen/Core>
#include <vector>

#include "gelex/exception.h"
#include "gelex/infra/record_visitor.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/state_record_set.h"

namespace gelex::bayes
{

struct ProportionState
{
    ProportionState() = default;
    ProportionState(
        const MixtureProportion& proportion,
        Eigen::Index num_markers);

    Eigen::VectorXi assignment;
    Eigen::VectorXi count;
    Eigen::VectorXd value;
    UpdatePolicy update{UpdatePolicy::fixed};
};

struct ComponentState
{
    ComponentState() = default;
    ComponentState(Eigen::Index num_components, Eigen::Index num_individuals);

    std::vector<Eigen::VectorXd> gebv;
    Eigen::VectorXd gebv_var;
};

class GeneticPriorState
{
   public:
    auto operator=(const GeneticPriorState&) -> GeneticPriorState& = delete;
    auto operator=(GeneticPriorState&&) noexcept -> GeneticPriorState& = delete;

    virtual ~GeneticPriorState() = default;

    template <typename Capability>
    auto query() -> Capability*
    {
        return dynamic_cast<Capability*>(this);
    }

    template <typename Capability>
    auto query() const -> const Capability*
    {
        return dynamic_cast<const Capability*>(this);
    }

    template <typename Capability>
    auto require() -> Capability&
    {
        auto* capability = query<Capability>();
        if (capability == nullptr)
        {
            throw GelexException(
                "genetic prior state lacks required capability");
        }
        return *capability;
    }

    template <typename Capability>
    auto require() const -> const Capability&
    {
        const auto* capability = query<Capability>();
        if (capability == nullptr)
        {
            throw GelexException(
                "genetic prior state lacks required capability");
        }
        return *capability;
    }

    virtual auto visit_records(StateRecordSet set, infra::RecordSink& sink)
        const -> void = 0;
    virtual auto visit_records(
        StateRecordSet set,
        infra::MutableRecordSink& sink) -> void = 0;

   protected:
    GeneticPriorState() = default;
    GeneticPriorState(const GeneticPriorState&) = default;
    GeneticPriorState(GeneticPriorState&&) noexcept = default;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_PRIOR_STATE_H_
