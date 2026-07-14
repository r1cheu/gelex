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

#ifndef GELEX_ALGO_REML_VARIANCE_CALCULATOR_H_
#define GELEX_ALGO_REML_VARIANCE_CALCULATOR_H_

#include <Eigen/Core>

#include "gelex/algo/reml/reml_buffer.h"
#include "gelex/freq/model.h"

namespace gelex
{

// Compute V = sum(K_i * sigma_i) + I * sigma_e
auto compute_v(
    const gelex::FreqModel& model,
    const gelex::FreqState& state,
    Eigen::Ref<Eigen::MatrixXd> v) -> void;

// In-place compute V^{-1} and return log|V| using Cholesky decomposition
// Input v is overwritten with V^{-1}
auto v_inv_logdet(Eigen::Ref<Eigen::MatrixXd> v) -> double;

// Compute the projection in factored form and Py = P*y.
// Fills ViX = V^{-1}*X, inv_XtViX = (X'*V^{-1}*X)^{-1}, logdet_xvx.
// P is represented lazily as V^{-1} - ViX*inv_XtViX*ViX' and never
// materialized. Requires v to already contain V^{-1}.
auto compute_proj(const gelex::FreqModel& model, RemlBuffer& buffer) -> void;

// Compute REML log-likelihood:
// logL = -0.5 * (log|V| + log|X'V^{-1}X| + y'Py)
auto compute_loglike(const gelex::FreqModel& model, const RemlBuffer& buffer)
    -> double;

// One O(n^3) evaluation at the variance components currently in `state`: fills
// V^{-1}, P, Py in buffer and returns the REML log-likelihood.
auto evaluate_point(
    const gelex::FreqModel& model,
    const gelex::FreqState& state,
    RemlBuffer& buffer) -> double;

// Pack variance components from FreqState into a vector (residual first).
auto collect_variance_components(const gelex::FreqState& state)
    -> Eigen::VectorXd;

// Unpack a variance-component vector back into FreqState (residual first).
auto distribute_variance_components(
    gelex::FreqState& state,
    const Eigen::Ref<const Eigen::VectorXd>& sigma) -> void;

}  // namespace gelex

#endif  // GELEX_ALGO_REML_VARIANCE_CALCULATOR_H_
