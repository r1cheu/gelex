#include "gelex/estimator/mcmc_result.h"

#include <algorithm>

#include "gelex/estimator/mcmc_samples.h"
#include "gelex/model/bayes.h"

namespace gelex
{
double compute_n_eff(const arma::dvec& samples)
{
    // Simplified effective sample size calculation
    size_t n = samples.n_elem;
    double mean_val = arma::mean(samples);
    double var_val = arma::var(samples);

    if (var_val == 0.0)
    {
        return static_cast<double>(n);
    }

    double sum_prod = 0.0;
    for (size_t i = 1; i < n; ++i)
    {
        sum_prod += (samples[i - 1] - mean_val) * (samples[i] - mean_val);
    }
    double autocorr = sum_prod / ((n - 1) * var_val);

    return static_cast<double>(n) / (1.0 + 2.0 * std::max(0.0, autocorr));
}

PosteriorStats compute_effect_result(const arma::dvec& samples)
{
    PosteriorStats result;
    result.means = {arma::mean(samples)};
    result.stds = {arma::stddev(samples)};
    result.medians = {arma::median(samples)};
    arma::dvec quantiles = arma::quantile(samples, arma::dvec{0.05, 0.95});
    result.q5s = {quantiles[0]};
    result.q95s = {quantiles[1]};
    result.n_effs = {compute_n_eff(samples)};
    result.r_hats = {1.0};  // Single chain
    return result;
}

PosteriorStats compute_effect_result(const arma::dmat& samples)
{
    PosteriorStats result;
    const size_t n_params = samples.n_rows;

    result.means.set_size(n_params);
    result.stds.set_size(n_params);
    result.medians.set_size(n_params);
    result.q5s.set_size(n_params);
    result.q95s.set_size(n_params);
    result.n_effs.set_size(n_params);
    result.r_hats.set_size(n_params);
    result.r_hats.fill(1.0);  // Single chain

    for (size_t i = 0; i < n_params; ++i)
    {
        arma::dvec param_samples = samples.row(i).t();
        result.means[i] = arma::mean(param_samples);
        result.stds[i] = arma::stddev(param_samples);
        result.medians[i] = arma::median(param_samples);
        arma::dvec quantiles
            = arma::quantile(param_samples, arma::dvec{0.05, 0.95});
        result.q5s[i] = quantiles[0];
        result.q95s[i] = quantiles[1];
        result.n_effs[i] = compute_n_eff(param_samples);
    }
    return result;
}

// MCMCResult compute_mcmc_result(
//     const MCMCSamples& samples,
//     const BayesModel& model)
// {
//     MCMCResult result;
//
//     // Mu effect (always scalar)
//     result.mu = compute_effect_result(samples.mu());
//
//     // Fixed effects
//     if (samples.fixed().n_rows > 0)
//     {
//         result.fixed = compute_effect_result(samples.fixed());
//     }
//
//     process_posterior(samples.random(), model.random(), result.random);
//     process_posterior(samples.genetic(), model.genetic(), result.genetic);
//
//     result.residual = compute_effect_result(samples.residual());
//
//     return result;
// }
}  // namespace gelex
