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

#include <string_view>

#include <Eigen/Core>

#include "gelex/model/bayes/capabilities.h"
#include "gelex/model/bayes/prior_state.h"

namespace gelex::bayes
{

class SingleMixtureAssignmentStateCap
{
   public:
    static constexpr std::string_view name = "single mixture assignment state";

    auto operator=(const SingleMixtureAssignmentStateCap&)
        -> SingleMixtureAssignmentStateCap& = delete;
    auto operator=(SingleMixtureAssignmentStateCap&&) noexcept
        -> SingleMixtureAssignmentStateCap& = delete;

    virtual ~SingleMixtureAssignmentStateCap() = default;

    virtual auto assignment() -> MixtureAssignmentState& = 0;
    virtual auto assignment() const -> const MixtureAssignmentState& = 0;

   protected:
    SingleMixtureAssignmentStateCap() = default;
    SingleMixtureAssignmentStateCap(const SingleMixtureAssignmentStateCap&)
        = default;
    SingleMixtureAssignmentStateCap(SingleMixtureAssignmentStateCap&&) noexcept
        = default;
};

class JointMixtureAssignmentStateCap
{
   public:
    static constexpr std::string_view name = "joint mixture assignment state";

    auto operator=(const JointMixtureAssignmentStateCap&)
        -> JointMixtureAssignmentStateCap& = delete;
    auto operator=(JointMixtureAssignmentStateCap&&) noexcept
        -> JointMixtureAssignmentStateCap& = delete;

    virtual ~JointMixtureAssignmentStateCap() = default;

    virtual auto assignment() -> MixtureAssignmentState& = 0;
    virtual auto assignment() const -> const MixtureAssignmentState& = 0;

   protected:
    JointMixtureAssignmentStateCap() = default;
    JointMixtureAssignmentStateCap(const JointMixtureAssignmentStateCap&)
        = default;
    JointMixtureAssignmentStateCap(JointMixtureAssignmentStateCap&&) noexcept
        = default;
};

class SingleComponentStateCap
{
   public:
    static constexpr std::string_view name = "single component state";

    auto operator=(const SingleComponentStateCap&)
        -> SingleComponentStateCap& = delete;
    auto operator=(SingleComponentStateCap&&) noexcept
        -> SingleComponentStateCap& = delete;

    virtual ~SingleComponentStateCap() = default;

    virtual auto component() -> ComponentState& = 0;
    virtual auto component() const -> const ComponentState& = 0;

   protected:
    SingleComponentStateCap() = default;
    SingleComponentStateCap(const SingleComponentStateCap&) = default;
    SingleComponentStateCap(SingleComponentStateCap&&) noexcept = default;
};

class JointComponentStateCap
{
   public:
    static constexpr std::string_view name = "joint component state";

    auto operator=(const JointComponentStateCap&)
        -> JointComponentStateCap& = delete;
    auto operator=(JointComponentStateCap&&) noexcept
        -> JointComponentStateCap& = delete;

    virtual ~JointComponentStateCap() = default;

    virtual auto component() -> ComponentState& = 0;
    virtual auto component() const -> const ComponentState& = 0;

   protected:
    JointComponentStateCap() = default;
    JointComponentStateCap(const JointComponentStateCap&) = default;
    JointComponentStateCap(JointComponentStateCap&&) noexcept = default;
};

using SingleSharedVarianceStateCap = SingleSharedVarianceCapability<double>;
using SinglePerMarkerVarianceStateCap
    = SinglePerMarkerVarianceCapability<Eigen::VectorXd>;
using SingleSampledMixtureProportionStateCap
    = SingleProportionCapability<SampledMixtureProportionState>;

using JointSharedVarianceStateCap = JointSharedVarianceCapability<double>;
using JointSampledMixtureProportionStateCap
    = JointProportionCapability<SampledMixtureProportionState>;

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_STATE_CAPABILITIES_H_
