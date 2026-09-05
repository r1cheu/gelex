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

#ifndef GELEX_BAYES_BASIC_STATE_H_
#define GELEX_BAYES_BASIC_STATE_H_

#include <Eigen/Core>

namespace gelex
{

struct FixedEffectState
{
    Eigen::VectorXd coefficients;
};

struct RandomEffectState
{
    Eigen::VectorXd coefficients;
    Eigen::VectorXd fitted_values;
    double variance{};
};

struct ResidualState
{
    Eigen::VectorXd adjusted_response;
    double variance{};
};

}  // namespace gelex

#endif  // GELEX_BAYES_BASIC_STATE_H_
