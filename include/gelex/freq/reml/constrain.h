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

#ifndef GELEX_FREQ_REML_CONSTRAIN_H_
#define GELEX_FREQ_REML_CONSTRAIN_H_

#include <Eigen/Dense>

namespace gelex
{
constexpr double constraint_scale = 1e-6;

// Lower limit that variance components are clamped to. A component at or below
// this floor is statistically indistinguishable from zero, so its Wald SE has
// no meaning.
inline auto constraint_floor(double y_variance) noexcept -> double
{
    return y_variance * constraint_scale;
}

// Clamps negative variance components to a small positive limit, redistributing
// the borrowed mass across the unconstrained ones. Returns a boundary mask
// flagging which components were clamped, so the caller can judge reliability
// (mask.count()) and mark boundary components whose Wald test is invalid.
auto constrain(Eigen::Ref<Eigen::VectorXd> varcmp, double y_variance)
    -> Eigen::ArrayX<bool>;
}  // namespace gelex

#endif  // GELEX_FREQ_REML_CONSTRAIN_H_
