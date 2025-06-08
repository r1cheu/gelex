#include "gelex/estimator/mcmc_samples.h"
#include "gelex/estimator/mcmc.h"
#include "gelex/model/bayes.h"
#include "gelex/model/bayes_effects.h"

namespace gelex
{

MCMCSamples::MCMCSamples(const MCMCParams& params, const BayesStatus& status, size_t n_chains)
    : n_records_((params.iter - params.n_burnin) / params.n_thin),
      n_chains_(n_chains)
{
    mu_.set_size(n_records_, n_chains_);
    if (status.fixed)
    {
        fixed_.set_size(status.fixed.coeff.n_elem, n_records_, n_chains_);
    }

    random_.coeffs.resize(status.random.size());
    random_.sigmas.resize(status.random.size());
    for (size_t i = 0; i < status.random.size(); ++i)
    {
        random_.coeffs[i].set_size(status.random[i].coeff.n_elem, n_records_, n_chains_);
        random_.sigmas[i].set_size(status.random[i].sigma.n_elem, n_records_, n_chains_);
    }

    genetic_.coeffs.resize(status.genetic.size());
    genetic_.sigmas.resize(status.genetic.size());
    for (size_t i = 0; i < status.genetic.size(); ++i)
    {
        genetic_.coeffs[i].set_size(status.genetic[i].coeff.n_elem, n_records_, n_chains_);
        genetic_.sigmas[i].set_size(status.genetic[i].sigma.n_elem, n_records_, n_chains_);
    }

    residual_.set_size(n_records_, n_chains_);
}

void MCMCSamples::store(const BayesStatus& status, size_t record_idx, size_t chain_idx)
{
    if (record_idx >= n_records_ || chain_idx >= n_chains_)
    {
        return;
    }

    mu_.at(record_idx, chain_idx) = status.mu.value;

    if (status.fixed)
    {
        fixed_.slice(chain_idx).col(record_idx) = status.fixed.coeff;
    }

    store_group(random_, status.random, record_idx, chain_idx);
    store_group(genetic_, status.genetic, record_idx, chain_idx);

    residual_.at(record_idx, chain_idx) = status.residual.value;
}

}  // namespace gelex
