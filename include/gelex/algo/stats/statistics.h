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

#ifndef GELEX_ESTIMATOR_FREQ_STATISTICS_H_
#define GELEX_ESTIMATOR_FREQ_STATISTICS_H_

namespace gelex
{

class FreqModel;
class FreqState;
class OptimizerState;

namespace statistics
{

// AIC = -2*logL + 2*k
// k = number of variance components + number of fixed effects
auto compute_aic(const FreqModel& model, double loglike) -> double;

// BIC = -2*logL + k*log(n)
auto compute_bic(const FreqModel& model, double loglike) -> double;

// Compute variance component standard errors from AI Hessian inverse
// se(σ) = sqrt(diag(-H⁻¹))
auto compute_variance_se(FreqState& state, const OptimizerState& opt_state)
    -> void;

// Compute heritability and its standard error using delta method
// h² = σ_g / Σσ
// se(h²) = sqrt(g' * (-H⁻¹) * g)
auto compute_heritability(FreqState& state, const OptimizerState& opt_state)
    -> void;

}  // namespace statistics
}  // namespace gelex

#endif  // GELEX_ESTIMATOR_FREQ_STATISTICS_H_
