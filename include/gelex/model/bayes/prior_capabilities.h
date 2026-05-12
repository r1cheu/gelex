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

#ifndef GELEX_MODEL_BAYES_PRIOR_CAPABILITIES_H_
#define GELEX_MODEL_BAYES_PRIOR_CAPABILITIES_H_

#include <span>

#include <Eigen/Core>

#include "gelex/model/bayes/prior_specs.h"

namespace gelex::bayes
{

class MarkerVarianceCapability
{
   public:
    auto operator=(const MarkerVarianceCapability&)
        -> MarkerVarianceCapability& = delete;
    auto operator=(MarkerVarianceCapability&&) noexcept
        -> MarkerVarianceCapability& = delete;

    virtual ~MarkerVarianceCapability() = default;

    virtual auto variance_specs() const -> std::span<const MarkerVarianceSpec>
        = 0;

   protected:
    MarkerVarianceCapability() = default;
    MarkerVarianceCapability(const MarkerVarianceCapability&) = default;
    MarkerVarianceCapability(MarkerVarianceCapability&&) noexcept = default;
};

class MixtureCapability
{
   public:
    auto operator=(const MixtureCapability&) -> MixtureCapability& = delete;
    auto operator=(MixtureCapability&&) noexcept -> MixtureCapability& = delete;

    virtual ~MixtureCapability() = default;

    virtual auto proportion_specs() const -> std::span<const ProportionSpec>
        = 0;

   protected:
    MixtureCapability() = default;
    MixtureCapability(const MixtureCapability&) = default;
    MixtureCapability(MixtureCapability&&) noexcept = default;
};

class ScaledMixtureCapability
{
   public:
    auto operator=(const ScaledMixtureCapability&)
        -> ScaledMixtureCapability& = delete;
    auto operator=(ScaledMixtureCapability&&) noexcept
        -> ScaledMixtureCapability& = delete;

    virtual ~ScaledMixtureCapability() = default;

    virtual auto multipliers() const -> std::span<const Eigen::VectorXd> = 0;

   protected:
    ScaledMixtureCapability() = default;
    ScaledMixtureCapability(const ScaledMixtureCapability&) = default;
    ScaledMixtureCapability(ScaledMixtureCapability&&) noexcept = default;
};

class JointMixtureCapability
{
   public:
    auto operator=(const JointMixtureCapability&)
        -> JointMixtureCapability& = delete;
    auto operator=(JointMixtureCapability&&) noexcept
        -> JointMixtureCapability& = delete;

    virtual ~JointMixtureCapability() = default;

    virtual auto proportion_spec() const -> const ProportionSpec& = 0;

   protected:
    JointMixtureCapability() = default;
    JointMixtureCapability(const JointMixtureCapability&) = default;
    JointMixtureCapability(JointMixtureCapability&&) noexcept = default;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_PRIOR_CAPABILITIES_H_
