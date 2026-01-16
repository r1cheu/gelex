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
