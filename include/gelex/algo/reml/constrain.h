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

#ifndef GELEX_ALGO_REML_CONSTRAIN_H_
#define GELEX_ALGO_REML_CONSTRAIN_H_

#include <Eigen/Dense>

namespace gelex
{
constexpr double CONSTRAINT_SCALE = 1e-6;

// Lower limit that variance components are clamped to. A component at or below
// this floor is statistically indistinguishable from zero, so its Wald SE has
// no meaning.
inline auto constraint_floor(double y_variance) noexcept -> double
{
    return y_variance * CONSTRAINT_SCALE;
}

// Clamps negative variance components to a small positive limit, redistributing
// the borrowed mass across the unconstrained ones. Returns how many components
// were constrained so the caller can judge the estimate's reliability.
auto constrain(Eigen::Ref<Eigen::VectorXd> varcmp, double y_variance)
    -> Eigen::Index;
}  // namespace gelex

#endif  // GELEX_ALGO_REML_CONSTRAIN_H_
