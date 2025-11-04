#include "gelex/types/mcmc_samples.h"

#include <ranges>

#include <Eigen/Core>

#include "gelex/estimator/bayes/params.h"
#include "gelex/model/bayes/model.h"
#include "types/bayes_effects.h"

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
RandomSamples::RandomSamples(
    const MCMCParams& params,
    const bayes::RandomEffect& effect)
    : RandomSamples(params, effect.design_matrix.cols()) {};

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

BaseMarkerSamples::BaseMarkerSamples(
    const MCMCParams& params,
    const bayes::GeneticEffect& effect)
    : RandomSamples(params, bayes::get_cols(effect.design_matrix))
{
    if (effect.init_pi)
    {
        const Eigen::Index num_snp = bayes::get_cols(effect.design_matrix);
        tracker.reserve(params.n_chains);
        n_props = effect.init_pi->size();
        for (Eigen::Index i = 0; i < params.n_chains; ++i)
        {
            tracker.emplace_back(num_snp, params.n_records);
        }
    }
    if (effect.estimate_pi)
    {
        const Eigen::Index n_props = effect.init_pi->size();
        prop.reserve(params.n_chains);
        for (Eigen::Index i = 0; i < params.n_chains; ++i)
        {
            prop.emplace_back(n_props, params.n_records);
        }
    }
}

AdditiveSamples::AdditiveSamples(
    const MCMCParams& params,
    const bayes::AdditiveEffect& effect)
    : BaseMarkerSamples(params, effect)
{
}

DominantSamples::DominantSamples(
    const MCMCParams& params,
    const bayes::DominantEffect& effect)
    : BaseMarkerSamples(params, effect)
{
}

ResidualSamples::ResidualSamples(const MCMCParams& params)
{
    variance.reserve(params.n_chains);
    for (Eigen::Index i = 0; i < params.n_chains; ++i)
    {
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

        if (!additive_->prop.empty())
        {
            additive_->prop[chain_idx].col(record_idx) = state->pi.prop;
        }
        if (!additive_->tracker.empty())
        {
            additive_->tracker[chain_idx].col(record_idx) = state->tracker;
        }
    }

    if (const auto* state = states.dominant(); dominant_ && state != nullptr)
    {
        dominant_->coeffs[chain_idx].col(record_idx) = state->coeffs;
        dominant_->variance[chain_idx](0, record_idx) = state->variance;

        if (!dominant_->prop.empty() && state->pi.prop.size() != 0)
        {
            dominant_->prop[chain_idx].col(record_idx) = state->pi.prop;
        }
        if (!dominant_->tracker.empty() && state->tracker.size() != 0)
        {
            dominant_->tracker[chain_idx].col(record_idx) = state->tracker;
        }
    }

    residual_.variance[chain_idx](0, record_idx) = states.residual().variance;
}

}  // namespace gelex
