#include "gelex/estimator/mcmc_samples.h"
#include "gelex/estimator/mcmc.h"
#include "gelex/model/bayes.h"
#include "gelex/model/bayes_effects.h"

namespace gelex
{

MCMCSamples::MCMCSamples(const MCMCParams& params, const BayesStatus& status)
    : n_records_((params.iter - params.n_burnin) / params.n_thin)
{
    mu_.set_size(n_records_);
    if (status.fixed)
    {
        fixed_.set_size(status.fixed.coeff.n_elem, n_records_);
    }

    init_group(random_, status.random);
    init_group(genetic_, status.genetic);

    residual_.set_size(n_records_);
}

void MCMCSamples::store(const BayesStatus& status, size_t record_idx)
{
    if (record_idx >= n_records_)
    {
        return;
    }

    mu_.at(record_idx) = status.mu.value;

    if (status.fixed)
    {
        fixed_.col(record_idx) = status.fixed.coeff;
    }

    store_group(random_, status.random, record_idx);
    store_group(genetic_, status.genetic, record_idx);

    residual_.at(record_idx) = status.residual.value;
}

}  // namespace gelex
