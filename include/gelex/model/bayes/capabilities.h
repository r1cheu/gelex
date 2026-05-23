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

#ifndef GELEX_MODEL_BAYES_CAPABILITIES_H_
#define GELEX_MODEL_BAYES_CAPABILITIES_H_

#include <span>

#include <Eigen/Core>

#include "gelex/model/bayes/prior_specs.h"

namespace gelex::bayes
{

template <typename T>
class VarianceCapability
{
   public:
    using element_type = T;

    auto operator=(const VarianceCapability&) -> VarianceCapability& = delete;
    auto operator=(VarianceCapability&&) noexcept
        -> VarianceCapability& = delete;

    virtual ~VarianceCapability() = default;

    virtual auto variance() -> std::span<T> = 0;
    virtual auto variance() const -> std::span<const T> = 0;

   protected:
    VarianceCapability() = default;
    VarianceCapability(const VarianceCapability&) = default;
    VarianceCapability(VarianceCapability&&) noexcept = default;
};

template <typename T>
class ProportionCapability
{
   public:
    using element_type = T;

    auto operator=(const ProportionCapability&)
        -> ProportionCapability& = delete;
    auto operator=(ProportionCapability&&) noexcept
        -> ProportionCapability& = delete;

    virtual ~ProportionCapability() = default;

    virtual auto proportion() -> std::span<T> = 0;
    virtual auto proportion() const -> std::span<const T> = 0;

   protected:
    ProportionCapability() = default;
    ProportionCapability(const ProportionCapability&) = default;
    ProportionCapability(ProportionCapability&&) noexcept = default;
};

template <typename T>
class MultiplierCapability
{
   public:
    using element_type = T;

    auto operator=(const MultiplierCapability&)
        -> MultiplierCapability& = delete;
    auto operator=(MultiplierCapability&&) noexcept
        -> MultiplierCapability& = delete;

    virtual ~MultiplierCapability() = default;

    virtual auto multiplier() -> std::span<T> = 0;
    virtual auto multiplier() const -> std::span<const T> = 0;

   protected:
    MultiplierCapability() = default;
    MultiplierCapability(const MultiplierCapability&) = default;
    MultiplierCapability(MultiplierCapability&&) noexcept = default;
};

using VarianceSpecCap = VarianceCapability<MarkerVarianceSpec>;
using ProportionSpecCap = ProportionCapability<ProportionSpec>;
using MultiplierSpecCap = MultiplierCapability<Eigen::VectorXd>;

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_CAPABILITIES_H_
