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

#ifndef GELEX_BAYES_GENETIC_PARAMETER_H_
#define GELEX_BAYES_GENETIC_PARAMETER_H_

#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>

#include "gelex/bayes/genetic_family.h"

namespace gelex
{

struct BetaHyperPrior
{
    double alpha{1.0};
    double beta{1.0};
};

template <std::size_t K>
struct DirichletHyperPrior
{
    constexpr DirichletHyperPrior() { concentration.fill(1.0); }

    explicit constexpr DirichletHyperPrior(std::array<double, K> concentration)
        : concentration{std::move(concentration)}
    {
    }

    std::array<double, K> concentration;
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

template <MixtureWeightUpdate Update>
using ProbabilityParameter = std::conditional_t<
    Update == MixtureWeightUpdate::Disabled,
    FixedParameter<double>,
    SampledParameter<double, BetaHyperPrior>>;

template <std::size_t ClassCount, MixtureWeightUpdate Update>
using SimplexParameter = std::conditional_t<
    Update == MixtureWeightUpdate::Disabled,
    FixedParameter<std::array<double, ClassCount>>,
    SampledParameter<
        std::array<double, ClassCount>,
        DirichletHyperPrior<ClassCount>>>;

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_PARAMETER_H_
