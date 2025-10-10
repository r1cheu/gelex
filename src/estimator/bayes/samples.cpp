#include "gelex/estimator/bayes/samples.h"

#include <Eigen/Core>

#include "../src/model/bayes/bayes_effects.h"
#include "gelex/estimator/bayes/params.h"
#include "gelex/model/bayes/model.h"

namespace gelex
{
using Eigen::Index;

FixedSamples FixedSamples::create(
    const MCMCParams& params,
    const BayesModel& model)
{
    FixedSamples samples;
    samples.coeffs.reserve(params.n_chains);
    samples.names = model.fixed()->names;
    const Eigen::Index n_coeffs = model.fixed()->design_matrix.cols();
    for (Eigen::Index i = 0; i < params.n_chains; ++i)
    {
        samples.coeffs.emplace_back(n_coeffs, params.n_records);
    }
    return samples;
}

AdditiveSamples AdditiveSamples::create(
    const MCMCParams& params,
    const BayesModel& model)
{
    AdditiveSamples samples;
    samples.reserve(params.n_chains);
    samples.name = model.additive()->name;

    const Eigen::Index n_coeffs = model.additive()->design_matrix.cols();
    const Eigen::Index n_variance = model.additive()->marker_variance.size();
    for (Eigen::Index i = 0; i < params.n_chains; ++i)
    {
        samples.coeffs.emplace_back(n_coeffs, params.n_records);
        samples.effect_variance.emplace_back(n_variance, params.n_records);
        samples.marker_variance.emplace_back(1, params.n_records);
    }
    return samples;
}

DominantSamples DominantSamples::create(
    const MCMCParams& params,
    const BayesModel& model)
{
    DominantSamples samples;
    samples.reserve(params.n_chains);

    const Eigen::Index n_coeffs = model.dominant()->design_matrix.cols();
    for (Eigen::Index i = 0; i < params.n_chains; ++i)
    {
        samples.coeffs.emplace_back(n_coeffs, params.n_records);
        samples.effect_variance.emplace_back(1, params.n_records);
    }
    return samples;
}

MultiRandomSamples MultiRandomSamples::create(
    const MCMCParams& params,
    const BayesModel& model)
{
    MultiRandomSamples samples;
    samples.reserve(model.random().size());
    for (const auto& effect : model.random())
    {
        RandomSamples effect_samples;
        effect_samples.reserve(params.n_chains);
        effect_samples.name = effect.name;
        const Eigen::Index n_coeff = effect.design_matrix.cols();
        for (Eigen::Index i = 0; i < params.n_chains; ++i)
        {
            effect_samples.coeffs.emplace_back(n_coeff, params.n_records);
            effect_samples.effect_variance.emplace_back(1, params.n_records);
        }
        samples.chain_samples.push_back(std::move(effect_samples));
    }
    return samples;
}

ResidualSamples ResidualSamples::create(
    const MCMCParams& params,
    const BayesModel& /* model */)
{
    ResidualSamples samples;
    samples.reserve(params.n_chains);
    for (Eigen::Index i = 0; i < params.n_chains; ++i)
    {
        samples.variance.emplace_back(1, params.n_records);
    }
    return samples;
}

PiSamples PiSamples::create(const MCMCParams& params, const BayesModel& model)
{
    PiSamples samples;
    samples.reserve(params.n_chains);
    const Eigen::Index n_props = model.additive()->pi.size();
    for (Eigen::Index i = 0; i < params.n_chains; ++i)
    {
        samples.prop.emplace_back(n_props, params.n_records);
    }
    return samples;
}

MCMCSamples::MCMCSamples(const MCMCParams& params, const BayesModel& model)
{
    if (model.fixed())
    {
        store_fixed_ = true;
        fixed_ = FixedSamples::create(params, model);
    }

    if (model.random().size() > 0)
    {
        store_random_ = true;
        random_ = MultiRandomSamples::create(params, model);
    }

    additive_ = AdditiveSamples::create(params, model);
    if (model.trait()->estimate_pi())
    {
        store_pi_ = true;
        pi_ = PiSamples::create(params, model);
    }
    if (model.dominant())
    {
        store_dominant_ = true;
        dominant_ = DominantSamples::create(params, model);
    }
    residual_ = ResidualSamples::create(params, model);
}

void MCMCSamples::store(
    const BayesStatus& status,
    Eigen::Index record_idx,
    Eigen::Index chain_idx)
{
    if (store_fixed_)
    {
        fixed_.coeffs[chain_idx].col(record_idx) = status.fixed->coeff;
    }

    if (store_random_)
    {
        for (Index effect_idx = 0; effect_idx < status.random.size();
             ++effect_idx)
        {
            random_.chain_samples[effect_idx].coeffs[chain_idx].col(record_idx)
                = status.random[effect_idx].coeff;
            random_.chain_samples[effect_idx].effect_variance[chain_idx].col(
                record_idx)
                = status.random[effect_idx].effect_variance;
        }
    }

    additive_.coeffs[chain_idx].col(record_idx) = status.additive->coeff;
    additive_.effect_variance[chain_idx](0, record_idx)
        = status.additive->effect_variance;
    additive_.marker_variance[chain_idx].col(record_idx)
        = status.additive->marker_variance;

    if (store_pi_)
    {
        pi_.prop[chain_idx].col(record_idx) = status.additive->pi.prop;
    }

    if (store_dominant_)
    {
        dominant_.coeffs[chain_idx].col(record_idx) = status.dominant->coeff;
        dominant_.effect_variance[chain_idx](0, record_idx)
            = status.dominant->effect_variance;
    }

    residual_.variance[chain_idx](0, record_idx) = status.residual.value;
}

}  // namespace gelex
