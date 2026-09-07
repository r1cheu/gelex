// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_FREQ_REML_STATISTICS_H_
#define GELEX_FREQ_REML_STATISTICS_H_

#include <cstddef>

#include "gelex/freq/reml/summary.h"

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

#endif  // GELEX_FREQ_REML_STATISTICS_H_
