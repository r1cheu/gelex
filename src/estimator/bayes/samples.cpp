#include "gelex/estimator/bayes/samples.h"

#include <Eigen/Core>

#include "../src/model/bayes/bayes_effects.h"
#include "gelex/estimator/bayes/params.h"
#include "gelex/model/bayes/model.h"

namespace gelex
{
using Eigen::Index;

MCMCSamples::MCMCSamples(const MCMCParams& params, const BayesModel& model)
    : n_records_((params.n_iters - params.n_burnin) / params.n_thin),
      n_chains_(params.n_chains),
      store_pi_(model.trait()->estimate_pi())
{
    // Pre-allocate all matrices upfront for better performance
    fixed_.coeff.reserve(n_chains_);
    fixed_.names = model.fixed().names;

    residual_.reserve(n_chains_);

    if (store_pi_)
    {
        pi_.reserve(n_chains_);
        const Index pi_size = model.additive()->pi.size();
        for (Index i = 0; i < n_chains_; i++)
        {
            pi_.emplace_back(pi_size, n_records_);
        }
    }

    if (model.additive())
    {
        const Index n_coeffs = model.additive()->design_matrix.cols();
        const Index sigma_size = model.additive()->sigma.size();

        additive_.coeffs.reserve(n_chains_);
        additive_.sigmas.reserve(n_chains_);
        additive_.variance.reserve(n_chains_);

        for (Index i = 0; i < n_chains_; i++)
        {
            additive_.coeffs.emplace_back(n_coeffs, n_records_);
            additive_.sigmas.emplace_back(sigma_size, n_records_);
            additive_.variance.emplace_back(1, n_records_);
        }
    }

    if (model.dominant())
    {
        const Index n_coeffs = model.dominant()->design_matrix.cols();
        dominant_.coeffs.reserve(n_chains_);
        dominant_.variance.reserve(n_chains_);

        for (Index i = 0; i < n_chains_; i++)
        {
            dominant_.coeffs.emplace_back(n_coeffs, n_records_);
            dominant_.variance.emplace_back(1, n_records_);
        }
    }

    // Pre-allocate fixed and residual effects
    const Index n_coeffs = model.fixed().design_matrix.cols();
    for (Index i = 0; i < n_chains_; ++i)
    {
        fixed_.coeff.emplace_back(n_coeffs, n_records_);
        residual_.emplace_back(1, n_records_);
    }
    init_random(model.random());
}

void MCMCSamples::store(
    const BayesStatus& status,
    Eigen::Index record_idx,
    Eigen::Index chain_idx)
{
    fixed_.coeff[chain_idx].col(record_idx) = status.fixed.coeff;

    if (status.additive)
    {
        additive_.coeffs[chain_idx].col(record_idx) = status.additive->coeff;
        additive_.sigmas[chain_idx].col(record_idx) = status.additive->sigma;
        additive_.variance[chain_idx](0, record_idx)
            = status.additive->variance;

        if (store_pi_)
        {
            pi_[chain_idx].col(record_idx) = status.additive->pi.prop;
        }
    }

    if (status.dominant)
    {
        dominant_.coeffs[chain_idx].col(record_idx) = status.dominant->coeff;
        dominant_.variance[chain_idx](0, record_idx)
            = status.dominant->variance;
    }

    residual_[chain_idx](0, record_idx) = status.residual.value;
}

void MCMCSamples::init_random(const bayes::RandomEffectManager& effects)
{
    const auto n_effects = static_cast<Index>(effects.size());
    if (n_effects == 0)
    {
        return;
    }

    random_.reserve(n_effects);

    for (Index effect_idx = 0; effect_idx < n_effects; ++effect_idx)
    {
        RandomSamples effect_samples;
        effect_samples.coeffs.reserve(n_chains_);
        effect_samples.sigmas.reserve(n_chains_);
        effect_samples.name = effects[effect_idx].name;

        const Index n_coeff = effects[effect_idx].design_matrix.cols();

        for (Index chain_idx = 0; chain_idx < n_chains_; ++chain_idx)
        {
            effect_samples.coeffs.emplace_back(n_coeff, n_records_);
            effect_samples.sigmas.emplace_back(1, n_records_);
        }

        random_.chain_samples.push_back(std::move(effect_samples));
    }
}

void MCMCSamples::store_random(
    const std::vector<bayes::RandomStatus>& status,
    Eigen::Index record_idx,
    Eigen::Index chain_idx)
{
    const auto n_effects = static_cast<Index>(status.size());
    for (Index effect_idx = 0; effect_idx < n_effects; ++effect_idx)
    {
        random_.chain_samples[effect_idx].coeffs[chain_idx].col(record_idx)
            = status[effect_idx].coeff;
        random_.chain_samples[effect_idx].sigmas[chain_idx].col(record_idx)
            = status[effect_idx].sigma;
    }
}

}  // namespace gelex
