#include "gelex/estimator/bayes/samples.h"

#include "armadillo"

#include "gelex/estimator/bayes/mcmc.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/policy.h"

namespace gelex
{

MCMCSamples::MCMCSamples(const MCMCParams& params, const BayesModel& model)
    : n_records_((params.n_iters - params.n_burnin) / params.n_thin),
      n_chains_(params.n_chains)
{
    fixed_.set_size(model.fixed().design_mat.n_cols, n_records_, n_chains_);
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

    fixed_.slice(chain_idx).col(record_idx) = status.fixed.coeff;
    store_group(status.random, record_idx, chain_idx);
    store_group(status.genetic, record_idx, chain_idx);
    residual_.at(0, record_idx, chain_idx) = status.residual.value;
}

void MCMCSamples::init_group(
    RandomGroup& group,
    const RandomEffectDesignManager& effects) const
{
    auto n_effects = effects.size();
    group.coeffs.resize(n_effects);
    group.sigmas.resize(n_effects);
    for (size_t i = 0; i < n_effects; ++i)
    {
        group.coeffs[i].set_size(
            effects[i].design_mat.n_cols, n_records_, n_chains_);
        group.sigmas[i].set_size(
            effects[i].sigma.n_elem, n_records_, n_chains_);
    }
}

void MCMCSamples::init_group(
    GeneticGroup& group,
    const GeneticEffectDesignManager& effects) const
{
    auto n_effects = effects.size();
    group.coeffs.resize(n_effects);
    group.sigmas.resize(n_effects);
    group.genetic_var.resize(n_effects);
    group.heritability.resize(n_effects);
    group.pi.resize(effects.size());
    for (size_t i = 0; i < effects.size(); ++i)
    {
        group.coeffs[i].set_size(
            effects[i].design_mat.n_cols, n_records_, n_chains_);
        group.sigmas[i].set_size(
            effects[i].sigma.n_elem, n_records_, n_chains_);
        group.genetic_var[i].set_size(1, n_records_, n_chains_);
        group.heritability[i].set_size(1, n_records_, n_chains_);
        if (bayes_trait_estimate_pi[to_index(effects[i].type)])
        {
            group.pi[i].set_size(effects[i].pi.n_elem, n_records_, n_chains_);
        }
    }
}

void MCMCSamples::store_group(
    const std::vector<RandomEffectState>& status,
    size_t record_idx,
    size_t chain_idx)
{
    for (size_t i = 0; i < status.size(); ++i)
    {
        random_.coeffs[i].slice(chain_idx).col(record_idx) = status[i].coeff;
        random_.sigmas[i].slice(chain_idx).col(record_idx) = status[i].sigma;
    }
}

void MCMCSamples::store_group(
    const std::vector<GeneticEffectState>& status,
    size_t record_idx,
    size_t chain_idx)
{
    for (size_t i = 0; i < status.size(); ++i)
    {
        genetic_.coeffs[i].slice(chain_idx).col(record_idx) = status[i].coeff;
        genetic_.sigmas[i].slice(chain_idx).col(record_idx) = status[i].sigma;
        genetic_.genetic_var[i].slice(chain_idx).col(record_idx)
            = status[i].genetic_var;
        genetic_.heritability[i].slice(chain_idx).col(record_idx)
            = status[i].heritability;
        if (bayes_trait_estimate_pi[to_index(status[i].type)])
        {
            genetic_.pi[i].slice(chain_idx).col(record_idx) = status[i].pi.prop;
        }
    }
}
}  // namespace gelex
