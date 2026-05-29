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
#include <string_view>

#include <Eigen/Core>

#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

template <typename T>
class ProportionCapability
{
   public:
    static constexpr std::string_view name = "proportion";
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
    static constexpr std::string_view name = "multiplier";
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

using FixedMixtureProportionCap = ProportionCapability<FixedMixtureProportion>;
using SampledMixtureProportionCap
    = ProportionCapability<SampledMixtureProportion>;
using MultiplierCap = MultiplierCapability<Eigen::VectorXd>;

template <typename T>
class SingleSharedVarCapability
{
   public:
    static constexpr std::string_view name = "single_shared_variance";
    using element_type = T;

    auto operator=(const SingleSharedVarCapability&)
        -> SingleSharedVarCapability& = delete;
    auto operator=(SingleSharedVarCapability&&) noexcept
        -> SingleSharedVarCapability& = delete;

    virtual ~SingleSharedVarCapability() = default;

    virtual auto variance() -> T& = 0;
    virtual auto variance() const -> const T& = 0;

   protected:
    SingleSharedVarCapability() = default;
    SingleSharedVarCapability(const SingleSharedVarCapability&) = default;
    SingleSharedVarCapability(SingleSharedVarCapability&&) noexcept = default;
};

template <typename T>
class SinglePerMarkerVarCapability
{
   public:
    static constexpr std::string_view name = "single_per_marker_variance";
    using element_type = T;

    auto operator=(const SinglePerMarkerVarCapability&)
        -> SinglePerMarkerVarCapability& = delete;
    auto operator=(SinglePerMarkerVarCapability&&) noexcept
        -> SinglePerMarkerVarCapability& = delete;

    virtual ~SinglePerMarkerVarCapability() = default;

    virtual auto variance() -> T& = 0;
    virtual auto variance() const -> const T& = 0;

   protected:
    SinglePerMarkerVarCapability() = default;
    SinglePerMarkerVarCapability(const SinglePerMarkerVarCapability&) = default;
    SinglePerMarkerVarCapability(SinglePerMarkerVarCapability&&) noexcept
        = default;
};

template <typename T>
class SingleProportionCapability
{
   public:
    static constexpr std::string_view name = "single_proportion";
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
    static constexpr std::string_view name = "single_multiplier";
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
class JointSharedVarCapability
{
   public:
    static constexpr std::string_view name = "joint_shared_variance";
    using element_type = T;

    auto operator=(const JointSharedVarCapability&)
        -> JointSharedVarCapability& = delete;
    auto operator=(JointSharedVarCapability&&) noexcept
        -> JointSharedVarCapability& = delete;

    virtual ~JointSharedVarCapability() = default;

    virtual auto variance(GeneticMode mode) -> T& = 0;
    virtual auto variance(GeneticMode mode) const -> const T& = 0;

   protected:
    JointSharedVarCapability() = default;
    JointSharedVarCapability(const JointSharedVarCapability&) = default;
    JointSharedVarCapability(JointSharedVarCapability&&) noexcept = default;
};

template <typename T>
class JointProportionCapability
{
   public:
    static constexpr std::string_view name = "joint_proportion";
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

using SingleSharedMarkerVarCap
    = SingleSharedVarCapability<SharedMarkerVariance>;
using SinglePerMarkerVarCap = SinglePerMarkerVarCapability<PerMarkerVariance>;
using SingleFixedMixtureProportionCap
    = SingleProportionCapability<FixedMixtureProportion>;
using SingleSampledMixtureProportionCap
    = SingleProportionCapability<SampledMixtureProportion>;
using SingleMultiplierCap = SingleMultiplierCapability<Eigen::VectorXd>;

using JointSharedMarkerVarCap = JointSharedVarCapability<SharedMarkerVariance>;
using JointFixedMixtureProportionCap
    = JointProportionCapability<FixedMixtureProportion>;
using JointSampledMixtureProportionCap
    = JointProportionCapability<SampledMixtureProportion>;

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_CAPABILITIES_H_
