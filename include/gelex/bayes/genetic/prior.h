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
#include <optional>
#include <utility>

#include "gelex/bayes/semantic_method.h"

namespace gelex
{

struct ScaledInvChiSqPrior
{
    double degrees_of_freedom{-2.0};
    double scale{0.0};
};

struct BetaPrior
{
    double alpha{1.0};
    double beta{1.0};
};

template <std::size_t Classes>
struct DirichletPrior
{
    constexpr DirichletPrior() { concentration.fill(1.0); }

    explicit constexpr DirichletPrior(std::array<double, Classes> concentration)
        : concentration{std::move(concentration)}
    {
    }

    std::array<double, Classes> concentration;
};

// A quantity the user may pin. An absent hyperprior is UpdatePolicy::Fixed
// compiled down: nothing samples the value, so there is no prior to hold.
template <typename Value, typename Hyper>
struct Updatable
{
    Value initial;
    std::optional<Hyper> prior;
};

// A quantity the chain always samples: the calibrated starting value and the
// hyperprior it is drawn under.
struct VarianceParameter
{
    double initial{};
    ScaledInvChiSqPrior prior;
};

using ProbabilityParameter = Updatable<double, BetaPrior>;

template <std::size_t Classes>
using SimplexParameter
    = Updatable<std::array<double, Classes>, DirichletPrior<Classes>>;

template <Variance Kind>
struct GaussianPrior
{
    VarianceParameter variance;
};

template <Variance Kind>
struct SpikeSlabPrior
{
    VarianceParameter variance;
    ProbabilityParameter probability;
};

struct ScaledMixturePrior
{
    VarianceParameter variance;
    SimplexParameter<5> probabilities;
    std::array<double, 5> scales{};
};

struct JointSpikeSlabPrior
{
    SimplexParameter<4> probabilities;
    ProbabilityParameter positive_probability;
};

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_PRIOR_H_
