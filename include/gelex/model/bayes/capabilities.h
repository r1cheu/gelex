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

#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

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

using MixtureProportionCap = ProportionCapability<MixtureProportion>;
using MultiplierCap = MultiplierCapability<Eigen::VectorXd>;

template <typename T>
class SharedVarianceCapability
{
   public:
    using element_type = T;

    auto operator=(const SharedVarianceCapability&)
        -> SharedVarianceCapability& = delete;
    auto operator=(SharedVarianceCapability&&) noexcept
        -> SharedVarianceCapability& = delete;

    virtual ~SharedVarianceCapability() = default;

    virtual auto variance() -> T& = 0;
    virtual auto variance() const -> const T& = 0;

   protected:
    SharedVarianceCapability() = default;
    SharedVarianceCapability(const SharedVarianceCapability&) = default;
    SharedVarianceCapability(SharedVarianceCapability&&) noexcept = default;
};

template <typename T>
class PerMarkerVarianceCapability
{
   public:
    using element_type = T;

    auto operator=(const PerMarkerVarianceCapability&)
        -> PerMarkerVarianceCapability& = delete;
    auto operator=(PerMarkerVarianceCapability&&) noexcept
        -> PerMarkerVarianceCapability& = delete;

    virtual ~PerMarkerVarianceCapability() = default;

    virtual auto variance() -> T& = 0;
    virtual auto variance() const -> const T& = 0;

   protected:
    PerMarkerVarianceCapability() = default;
    PerMarkerVarianceCapability(const PerMarkerVarianceCapability&) = default;
    PerMarkerVarianceCapability(PerMarkerVarianceCapability&&) noexcept
        = default;
};

template <typename T>
class SingleProportionCapability
{
   public:
    using element_type = T;

    auto operator=(const SingleProportionCapability&)
        -> SingleProportionCapability& = delete;
    auto operator=(SingleProportionCapability&&) noexcept
        -> SingleProportionCapability& = delete;

    virtual ~SingleProportionCapability() = default;

    virtual auto proportion() -> T& = 0;
    virtual auto proportion() const -> const T& = 0;

   protected:
    SingleProportionCapability() = default;
    SingleProportionCapability(const SingleProportionCapability&) = default;
    SingleProportionCapability(SingleProportionCapability&&) noexcept = default;
};

template <typename T>
class SingleMultiplierCapability
{
   public:
    using element_type = T;

    auto operator=(const SingleMultiplierCapability&)
        -> SingleMultiplierCapability& = delete;
    auto operator=(SingleMultiplierCapability&&) noexcept
        -> SingleMultiplierCapability& = delete;

    virtual ~SingleMultiplierCapability() = default;

    virtual auto multiplier() -> T& = 0;
    virtual auto multiplier() const -> const T& = 0;

   protected:
    SingleMultiplierCapability() = default;
    SingleMultiplierCapability(const SingleMultiplierCapability&) = default;
    SingleMultiplierCapability(SingleMultiplierCapability&&) noexcept = default;
};

template <typename T>
class JointVarianceCapability
{
   public:
    using element_type = T;

    auto operator=(const JointVarianceCapability&)
        -> JointVarianceCapability& = delete;
    auto operator=(JointVarianceCapability&&) noexcept
        -> JointVarianceCapability& = delete;

    virtual ~JointVarianceCapability() = default;

    virtual auto variance(GeneticMode mode) -> T& = 0;
    virtual auto variance(GeneticMode mode) const -> const T& = 0;

   protected:
    JointVarianceCapability() = default;
    JointVarianceCapability(const JointVarianceCapability&) = default;
    JointVarianceCapability(JointVarianceCapability&&) noexcept = default;
};

template <typename T>
class JointProportionCapability
{
   public:
    using element_type = T;

    auto operator=(const JointProportionCapability&)
        -> JointProportionCapability& = delete;
    auto operator=(JointProportionCapability&&) noexcept
        -> JointProportionCapability& = delete;

    virtual ~JointProportionCapability() = default;

    virtual auto proportion() -> T& = 0;
    virtual auto proportion() const -> const T& = 0;

   protected:
    JointProportionCapability() = default;
    JointProportionCapability(const JointProportionCapability&) = default;
    JointProportionCapability(JointProportionCapability&&) noexcept = default;
};

using SingleSharedMarkerVarianceCap
    = SharedVarianceCapability<SharedMarkerVariance>;
using SinglePerMarkerVarianceCap
    = PerMarkerVarianceCapability<PerMarkerVariance>;
using SingleMixtureProportionCap
    = SingleProportionCapability<MixtureProportion>;
using SingleMultiplierCap = SingleMultiplierCapability<Eigen::VectorXd>;
using JointMarkerVarianceCap = JointVarianceCapability<JointMarkerVariance>;
using JointMixtureProportionCap = JointProportionCapability<MixtureProportion>;

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_CAPABILITIES_H_
