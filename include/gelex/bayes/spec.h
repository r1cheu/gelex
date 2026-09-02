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

namespace gelex
{

struct NoParameters
{
};

struct SpikeSlab
{
    double probability{0.01};
};

struct ScaledMixture
{
    std::array<double, 5> probabilities{0.99, 0.005, 0.003, 0.001, 0.001};
    std::array<double, 5> scales{0.0, 0.001, 0.01, 0.1, 1.0};
};

struct JointSpikeSlab
{
    std::array<double, 4> probabilities{0.99, 1.0 / 300, 1.0 / 300, 1.0 / 300};
    double positive_probability{0.5};
};

}  // namespace gelex

#endif  // GELEX_BAYES_SPEC_H_
