#include "gelex/estimator/bayes/samples.h"

#include <ranges>

#include <Eigen/Core>

#include "../src/model/bayes/bayes_effects.h"
#include "gelex/estimator/bayes/params.h"
#include "gelex/model/bayes/model.h"

namespace gelex
{
using Eigen::Index;

FixedSamples::FixedSamples(
    const MCMCParams& params,
    const bayes::FixedEffect& effect)
{
    coeffs.reserve(params.n_chains);
    for (Eigen::Index i = 0; i < params.n_chains; ++i)
    {
        coeffs.emplace_back(effect.design_matrix.cols(), params.n_records);
    }
}

RandomSamples::RandomSamples(const MCMCParams& params, Eigen::Index n_coeffs)
{
    coeffs.reserve(params.n_chains);
    variance.reserve(params.n_chains);
    for (Eigen::Index i = 0; i < params.n_chains; ++i)
    {
        coeffs.emplace_back(n_coeffs, params.n_records);
        variance.emplace_back(1, params.n_records);
    }
}

MCMCSamples::MCMCSamples(const MCMCParams& params, const BayesModel& model)
    : residual_(params)
{
    if (const auto* effect = model.fixed(); effect)
    {
        fixed_.emplace(params, *effect);
    }

    if (const auto& effects = model.random(); !effects.empty())
    {
        random_.reserve(effects.size());
        for (const auto& effect : effects)
        {
            random_.emplace_back(params, effect);
        }
    }

    if (const auto* effect = model.additive(); effect)
    {
        additive_.emplace(params, *effect);

        // Only allocate pi and tracker for mixture models
        if (model.is_mixture_model())
        {
            pi_.emplace(params, *effect);
            snp_tracker_.emplace(params, *effect);
        }
    }

    if (const auto* effect = model.dominant(); effect)
    {
        dominant_.emplace(params, *effect);
    }
}

void MCMCSamples::store(
    const BayesState& states,
    Eigen::Index record_idx,
    Eigen::Index chain_idx)
{
    if (const auto* state = states.fixed(); fixed_ && state != nullptr)
    {
        fixed_->coeffs[chain_idx].col(record_idx) = state->coeffs;
    }

    for (auto&& [sample, state] : std::views::zip(random_, states.random()))
    {
        sample.coeffs[chain_idx].col(record_idx) = state.coeffs;
        sample.variance[chain_idx](0, record_idx) = state.variance;
    }

    if (const auto* state = states.additive(); additive_ && state != nullptr)
    {
        additive_->coeffs[chain_idx].col(record_idx) = state->coeffs;
        additive_->variance[chain_idx](0, record_idx) = state->variance;

        // Store pi and tracker for mixture models
        if (pi_)
        {
            pi_->prop[chain_idx].col(record_idx) = state->pi.prop;
        }
        if (snp_tracker_)
        {
            snp_tracker_->tracker[chain_idx].col(record_idx) = state->tracker;
        }
    }

    if (const auto* state = states.dominant(); dominant_ && state != nullptr)
    {
        dominant_->coeffs[chain_idx].col(record_idx) = state->coeffs;
        dominant_->ratios[chain_idx].col(record_idx) = state->ratios;
        dominant_->variance[chain_idx](0, record_idx) = state->variance;
    }

    residual_.variance[chain_idx](0, record_idx) = states.residual().variance;
}

}  // namespace gelex
