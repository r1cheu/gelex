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

#ifndef GELEX_MODEL_BAYES_STATE_CAPABILITIES_H_
#define GELEX_MODEL_BAYES_STATE_CAPABILITIES_H_

#include <span>

#include <Eigen/Core>

#include "gelex/infra/record_visitor.h"
#include "gelex/model/bayes/capabilities.h"
#include "gelex/model/bayes/prior_state.h"

namespace gelex::bayes
{

class VarianceStateCap : public VarianceCapability<Eigen::VectorXd>
{
   protected:
    auto visit_sample_records(infra::RecordVisitor& visitor) const -> void;
    auto visit_checkpoint_records(infra::RecordVisitor& visitor) const -> void;
    auto visit_checkpoint_records(infra::MutableRecordVisitor& visitor) -> void;
};

class ProportionStateCap : public ProportionCapability<ProportionState>
{
   protected:
    auto visit_sample_records(infra::RecordVisitor& visitor) const -> void;
    auto visit_checkpoint_records(infra::RecordVisitor& visitor) const -> void;
    auto visit_checkpoint_records(infra::MutableRecordVisitor& visitor) -> void;
};

class ComponentStateCap
{
   public:
    auto operator=(const ComponentStateCap&) -> ComponentStateCap& = delete;
    auto operator=(ComponentStateCap&&) noexcept -> ComponentStateCap& = delete;

    virtual ~ComponentStateCap() = default;

    virtual auto component() -> std::span<ComponentState> = 0;
    virtual auto component() const -> std::span<const ComponentState> = 0;

   protected:
    ComponentStateCap() = default;
    ComponentStateCap(const ComponentStateCap&) = default;
    ComponentStateCap(ComponentStateCap&&) noexcept = default;

    auto visit_sample_records(infra::RecordVisitor& visitor) const -> void;
    auto visit_checkpoint_records(infra::RecordVisitor& visitor) const -> void;
    auto visit_checkpoint_records(infra::MutableRecordVisitor& visitor) -> void;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_STATE_CAPABILITIES_H_
