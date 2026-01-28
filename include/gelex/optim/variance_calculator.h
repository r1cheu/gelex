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

#ifndef GELEX_OPTIM_VARIANCE_CALCULATOR_H_
#define GELEX_OPTIM_VARIANCE_CALCULATOR_H_

#include <Eigen/Core>

namespace gelex
{

class FreqModel;
class FreqState;
class OptimizerState;

namespace variance_calculator
{

// Compute V = sum(Z_i * K_i * Z_i' * sigma_i) + I * sigma_e
// For genetic effects without Z, V += K * sigma
auto compute_v(
    const FreqModel& model,
    const FreqState& state,
    Eigen::Ref<Eigen::MatrixXd> v) -> void;

// In-place compute V^{-1} and return log|V| using Cholesky decomposition
// Input v is overwritten with V^{-1}
auto v_inv_logdet(Eigen::Ref<Eigen::MatrixXd> v) -> double;

// Compute projection matrix P = V^{-1} - V^{-1}*X*(X'*V^{-1}*X)^{-1}*X'*V^{-1}
// and Py = P * y
// Requires v to already contain V^{-1}
auto compute_proj(const FreqModel& model, OptimizerState& state) -> void;

// Compute REML log-likelihood:
// logL = -0.5 * (log|V| + log|X'V^{-1}X| + y'Py)
auto compute_loglike(const FreqModel& model, const OptimizerState& state)
    -> double;

}  // namespace variance_calculator
}  // namespace gelex

#endif  // GELEX_OPTIM_VARIANCE_CALCULATOR_H_
