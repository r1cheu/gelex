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

#include "gelex/bayes/genetic/parameter.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/variance/parameter.h"

namespace gelex
{

template <VarianceLayout Kind>
struct GaussianPrior
{
    VarianceParameter variance;
};

// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
struct HalfNormalPrior
{
    VarianceParameter variance;
    ProbabilityParameter<MixtureWeightUpdate::Enabled> positive_probability;
};

template <
    VarianceLayout Kind,
    MixtureWeightUpdate WeightUpdate = MixtureWeightUpdate::Enabled>
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
struct SpikeSlabPrior
{
    VarianceParameter variance;
    ProbabilityParameter<WeightUpdate> probability;
};

template <
    std::size_t ClassCount,
    MixtureWeightUpdate WeightUpdate = MixtureWeightUpdate::Enabled>
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
struct ScaledMixturePrior
{
    static constexpr std::size_t class_count = ClassCount;

    VarianceParameter variance;
    SimplexParameter<class_count, WeightUpdate> probabilities;
    std::array<double, class_count> scales{};
};

template <
    std::size_t ClassCount,
    MixtureWeightUpdate WeightUpdate = MixtureWeightUpdate::Enabled>
struct JointSpikeSlabPrior
{
    static constexpr std::size_t class_count = ClassCount;

    SimplexParameter<class_count, WeightUpdate> probabilities;
};

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_PRIOR_H_
