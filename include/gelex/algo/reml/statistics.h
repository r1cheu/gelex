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

#ifndef GELEX_ALGO_REML_STATISTICS_H_
#define GELEX_ALGO_REML_STATISTICS_H_

#include <cstddef>

#include "gelex/algo/reml/summary.h"

namespace gelex
{
class FreqModel;
class FreqState;
}  // namespace gelex

namespace gelex
{
class RemlBuffer;

// AIC = -2*logL + 2*k
// k = number of variance components + number of fixed effects
auto compute_aic(const FreqModel& model, double loglike) -> double;

// BIC = -2*logL + k*log(n)
auto compute_bic(const FreqModel& model, double loglike) -> double;

// One-sided Wald p-value P(Z > z) for a boundary parameter (variance
// components: H0 sigma^2 = 0 against sigma^2 > 0). NaN propagates unchanged.
auto wald_p_onesided(double z) noexcept -> double;

// Two-sided Wald p-value 2*P(Z > |z|) for an interior parameter (fixed
// effects). NaN propagates unchanged.
auto wald_p_twosided(double z) noexcept -> double;

// Compute variance component standard errors from AI Hessian inverse
// se(σ) = sqrt(diag(-H⁻¹))
auto compute_variance_se(FreqState& state, const RemlBuffer& buffer) -> void;

// Compute variance ratio and its standard error using delta method
// ratio = σ_r / Σσ
// se(ratio) = sqrt(g' * (-H⁻¹) * g)
auto compute_variance_ratio(FreqState& state, const RemlBuffer& buffer) -> void;

// Post-estimation: solve fixed/random effects, fill SEs and ratios, pack the
// reportable summary, and materialize P into buffer for the GWAS operators.
auto assemble_reml_fit(
    const FreqModel& model,
    FreqState& state,
    RemlBuffer& buffer,
    double loglike,
    bool converged,
    size_t iter_count) -> RemlFit;

}  // namespace gelex

#endif  // GELEX_ALGO_REML_STATISTICS_H_
