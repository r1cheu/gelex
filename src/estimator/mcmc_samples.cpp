#include "gelex/estimator/mcmc_samples.h"
#include "gelex/estimator/mcmc.h"
#include "gelex/model/bayes.h"
#include "gelex/model/bayes_effects.h"

namespace gelex
{

MCMCSamples::MCMCSamples(const MCMCParams& params, const BayesModel& model)
    : n_records_((params.iter - params.n_burnin) / params.n_thin),
      n_chains_(params.n_chains)
{
    mu_.set_size(1, n_records_, n_chains_);
    if (model.fixed())
    {
        fixed_.set_size(
            model.fixed()->design_mat.n_cols, n_records_, n_chains_);
    }
    init_group(random_, model.random());
    init_group(genetic_, model.genetic());
    residual_.set_size(1, n_records_, n_chains_);
}

void MCMCSamples::store(
    const BayesStatus& status,
    size_t record_idx,
    size_t chain_idx)
{
    if (record_idx >= n_records_ || chain_idx >= n_chains_)
    {
        return;
    }

    mu_.at(0, record_idx, chain_idx) = status.mu.value;

    if (status.fixed)
    {
        fixed_.slice(chain_idx).col(record_idx) = status.fixed.coeff;
    }

    store_group(random_, status.random, record_idx, chain_idx);
    store_group(genetic_, status.genetic, record_idx, chain_idx);

    residual_.at(0, record_idx, chain_idx) = status.residual.value;
}

}  // namespace gelex
