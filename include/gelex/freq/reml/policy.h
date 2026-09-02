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

#ifndef GELEX_FREQ_REML_POLICY_H_
#define GELEX_FREQ_REML_POLICY_H_

#include <Eigen/Core>

#include "gelex/freq/model.h"
#include "gelex/freq/reml/reml_buffer.h"

namespace gelex
{

struct EMPolicy
{
    static auto apply(
        const gelex::FreqModel& model,
        const gelex::FreqState& state,
        RemlBuffer& buffer) -> Eigen::VectorXd;
};

struct AIPolicy
{
    // AI-REML search direction delta = -H^{-1} * grad, evaluated at the point
    // whose V^{-1}/P/Py are currently held in buffer. The caller owns the
    // step length (see the Armijo backtracking in Estimator::fit).
    static auto direction(const gelex::FreqModel& model, RemlBuffer& buffer)
        -> Eigen::VectorXd;
};

}  // namespace gelex

#endif  // GELEX_FREQ_REML_POLICY_H_
