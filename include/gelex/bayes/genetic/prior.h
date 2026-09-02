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

#ifndef GELEX_BAYES_GENETIC_PRIOR_H_
#define GELEX_BAYES_GENETIC_PRIOR_H_

#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>

#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/variance_parameter.h"

namespace gelex
{

struct BetaHyperPrior
{
    double alpha{1.0};
    double beta{1.0};
};

template <std::size_t Classes>
struct DirichletHyperPrior
{
    constexpr DirichletHyperPrior() { concentration.fill(1.0); }

    explicit constexpr DirichletHyperPrior(
        std::array<double, Classes> concentration)
        : concentration{std::move(concentration)}
    {
    }

    std::array<double, Classes> concentration;
};

template <typename T>
struct FixedParameter
{
    T initial;
};

template <typename T, typename HyperPrior>
struct SampledParameter
{
    T initial;
    HyperPrior hyperprior;
};

template <UpdatePolicy Policy>
using ProbabilityParameter = std::conditional_t<
    Policy == UpdatePolicy::Fixed,
    FixedParameter<double>,
    SampledParameter<double, BetaHyperPrior>>;

template <std::size_t ClassCount, UpdatePolicy Policy>
using SimplexParameter = std::conditional_t<
    Policy == UpdatePolicy::Fixed,
    FixedParameter<std::array<double, ClassCount>>,
    SampledParameter<
        std::array<double, ClassCount>,
        DirichletHyperPrior<ClassCount>>>;

template <VarianceLayout Kind>
struct GaussianPrior
{
    VarianceParameter variance;
};

template <HalfNormalAsymmetry Asymmetry>
struct HalfNormalPrior;

template <>
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
struct HalfNormalPrior<HalfNormalAsymmetry::Count>
{
    VarianceParameter variance;
    ProbabilityParameter<UpdatePolicy::Sampled> positive_probability;
};

template <>
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
struct HalfNormalPrior<HalfNormalAsymmetry::Magnitude>
{
    VarianceParameter variance;
};

template <
    VarianceLayout Kind,
    UpdatePolicy ProbabilityUpdate = UpdatePolicy::Sampled>
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
struct SpikeSlabPrior
{
    VarianceParameter variance;
    ProbabilityParameter<ProbabilityUpdate> probability;
};

template <
    std::size_t ClassCount,
    UpdatePolicy ProbabilitiesUpdate = UpdatePolicy::Sampled>
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
struct ScaledMixturePrior
{
    static constexpr std::size_t class_count = ClassCount;

    VarianceParameter variance;
    SimplexParameter<class_count, ProbabilitiesUpdate> probabilities;
    std::array<double, class_count> scales{};
};

template <
    std::size_t ClassCount,
    UpdatePolicy ProbabilitiesUpdate = UpdatePolicy::Sampled>
struct JointSpikeSlabPrior
{
    static constexpr std::size_t class_count = ClassCount;

    SimplexParameter<class_count, ProbabilitiesUpdate> probabilities;
};

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_PRIOR_H_
