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

#ifndef GELEX_BAYES_SPEC_H_
#define GELEX_BAYES_SPEC_H_

#include <array>
#include <cstdint>

namespace gelex
{

enum class UpdatePolicy : std::uint8_t
{
    Fixed,
    Sampled,
};

template <typename Value>
struct MaybeSampled
{
    Value initial;
    UpdatePolicy update{UpdatePolicy::Sampled};
};

// The parameter shape of a method whose prior structure is fixed by the data
// alone. It occupies no storage and carries no value, so it states an absence
// at the type level rather than standing in for one at run time.
struct NoParameters
{
};

struct SpikeSlab
{
    MaybeSampled<double> probability{.initial = 0.01};
};

struct ScaledMixture
{
    MaybeSampled<std::array<double, 5>> probabilities{
        .initial = {0.99, 0.005, 0.003, 0.001, 0.001},
    };
    std::array<double, 5> scales{0.0, 0.001, 0.01, 0.1, 1.0};
};

struct JointSpikeSlab
{
    MaybeSampled<std::array<double, 4>> probabilities{
        .initial = {0.99, 1.0 / 300, 1.0 / 300, 1.0 / 300},
    };
    MaybeSampled<double> positive_probability{.initial = 0.5};
};

}  // namespace gelex

#endif  // GELEX_BAYES_SPEC_H_
