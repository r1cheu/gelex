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
