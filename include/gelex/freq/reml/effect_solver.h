// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_FREQ_REML_EFFECT_SOLVER_H_
#define GELEX_FREQ_REML_EFFECT_SOLVER_H_

namespace gelex
{
class FreqModel;
class FreqState;
}  // namespace gelex

namespace gelex
{
class RemlBuffer;

// Compute fixed effects coefficients (BLUE) and standard errors
// β = (X'V⁻¹X)⁻¹ X'V⁻¹y
// se(β) = sqrt(diag((X'V⁻¹X)⁻¹))
auto compute_fixed_effects(
    const FreqModel& model,
    FreqState& state,
    const RemlBuffer& buffer) -> void;

// Compute random effects (BLUP)
// u = K * Py * σ
auto compute_random_effects(
    const FreqModel& model,
    FreqState& state,
    const RemlBuffer& buffer) -> void;

}  // namespace gelex

#endif  // GELEX_FREQ_REML_EFFECT_SOLVER_H_
