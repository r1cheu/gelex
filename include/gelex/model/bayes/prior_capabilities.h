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

class VarianceCapability
{
   public:
    auto operator=(const VarianceCapability&) -> VarianceCapability& = delete;
    auto operator=(VarianceCapability&&) noexcept
        -> VarianceCapability& = delete;

    virtual ~VarianceCapability() = default;

    virtual auto variance_specs() const -> std::span<const MarkerVarianceSpec>
        = 0;

   protected:
    VarianceCapability() = default;
    VarianceCapability(const VarianceCapability&) = default;
    VarianceCapability(VarianceCapability&&) noexcept = default;
};

class ProportionCapability
{
   public:
    auto operator=(const ProportionCapability&)
        -> ProportionCapability& = delete;
    auto operator=(ProportionCapability&&) noexcept
        -> ProportionCapability& = delete;

    virtual ~ProportionCapability() = default;

    virtual auto proportion_specs() const -> std::span<const ProportionSpec>
        = 0;

   protected:
    ProportionCapability() = default;
    ProportionCapability(const ProportionCapability&) = default;
    ProportionCapability(ProportionCapability&&) noexcept = default;
};

class MultiplierCapability
{
   public:
    auto operator=(const MultiplierCapability&)
        -> MultiplierCapability& = delete;
    auto operator=(MultiplierCapability&&) noexcept
        -> MultiplierCapability& = delete;

    virtual ~MultiplierCapability() = default;

    virtual auto multipliers() const -> std::span<const Eigen::VectorXd> = 0;

   protected:
    MultiplierCapability() = default;
    MultiplierCapability(const MultiplierCapability&) = default;
    MultiplierCapability(MultiplierCapability&&) noexcept = default;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_PRIOR_CAPABILITIES_H_
