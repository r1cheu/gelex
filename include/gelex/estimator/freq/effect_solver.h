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

#ifndef GELEX_ESTIMATOR_FREQ_EFFECT_SOLVER_H_
#define GELEX_ESTIMATOR_FREQ_EFFECT_SOLVER_H_

namespace gelex
{

class FreqModel;
class FreqState;
class OptimizerState;

namespace effect_solver
{

// Compute fixed effects coefficients (BLUE) and standard errors
// β = (X'V⁻¹X)⁻¹ X'V⁻¹y
// se(β) = sqrt(diag((X'V⁻¹X)⁻¹))
auto compute_fixed_effects(
    const FreqModel& model,
    FreqState& state,
    const OptimizerState& opt_state) -> void;

// Compute random effects (BLUP)
// u = K * Py * σ
auto compute_random_effects(
    const FreqModel& model,
    FreqState& state,
    const OptimizerState& opt_state) -> void;

}  // namespace effect_solver
}  // namespace gelex

#endif  // GELEX_ESTIMATOR_FREQ_EFFECT_SOLVER_H_
