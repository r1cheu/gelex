#include "gelex/estimator/mcmc_result.h"

#include <algorithm>

#include "gelex/estimator/mcmc_storage.h"
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

    // Compute autocorrelation at lag 1
    double sum_prod = 0.0;
    for (size_t i = 1; i < n; ++i)
    {
        sum_prod += (samples[i - 1] - mean_val) * (samples[i] - mean_val);
    }
    double autocorr = sum_prod / ((n - 1) * var_val);

    // Simple effective sample size estimate
    return static_cast<double>(n) / (1.0 + 2.0 * std::max(0.0, autocorr));
}

ParameterResult compute_parameter_result(const arma::dvec& samples)
{
    ParameterResult result{};

    result.mean = arma::mean(samples);
    result.std = arma::stddev(samples);
    result.median = arma::median(samples);

    arma::dvec quantiles = arma::quantile(samples, arma::dvec{0.05, 0.95});
    result.q5 = quantiles[0];
    result.q95 = quantiles[1];

    result.n_eff = compute_n_eff(samples);
    result.r_hat = 1.0;  // Single chain, cannot compute R-hat

    return result;
}

EffectResult compute_effect_result(const arma::dmat& samples)
{
    EffectResult result;
    size_t n_params = samples.n_rows;
    result.parameters.resize(n_params);

    for (size_t i = 0; i < n_params; ++i)
    {
        arma::dvec param_samples = samples.row(i).t();
        result.parameters[i] = compute_parameter_result(param_samples);
    }

    return result;
}

EffectResult compute_effect_result(const arma::dvec& samples)
{
    EffectResult result;
    result.parameters.resize(1);
    result.parameters[0] = compute_parameter_result(samples);
    return result;
}

MCMCResult compute_mcmc_result(const MCMCStorage& storage, const Bayes& model)
{
    MCMCResult result;

    // Compute mu results
    result.mu = compute_effect_result(storage.mu_samples());

    // Compute fixed effect results
    if (model.fixed().exist)
    {
        result.fixed = compute_effect_result(storage.fixed_samples());
    }

    // Compute random effect results
    result.random.resize(storage.random_samples().size());
    result.random_sigma.resize(storage.random_sigma_samples().size());
    result.random_names.resize(model.random().size());

    for (size_t i = 0; i < storage.random_samples().size(); ++i)
    {
        result.random[i] = compute_effect_result(storage.random_samples()[i]);
        result.random_sigma[i]
            = compute_effect_result(storage.random_sigma_samples()[i]);
        result.random_names[i] = model.random()[i].name;
    }

    // Compute genetic effect results
    result.genetic.resize(storage.genetic_samples().size());
    result.genetic_sigma.resize(storage.genetic_sigma_samples().size());
    result.genetic_names.resize(model.genetic().size());

    for (size_t i = 0; i < storage.genetic_samples().size(); ++i)
    {
        result.genetic[i] = compute_effect_result(storage.genetic_samples()[i]);
        result.genetic_sigma[i]
            = compute_effect_result(storage.genetic_sigma_samples()[i]);
        result.genetic_names[i] = model.genetic()[i].name;
    }

    // Compute residual results
    result.residual = compute_effect_result(storage.residual_samples());

    return result;
}
}  // namespace gelex
