#include "gelex/estimator/mcmc_storage.h"
#include "gelex/estimator/mcmc.h"
#include "gelex/model/bayes.h"
#include "gelex/model/effects/bayes_effects.h"

namespace gelex
{

MCMCStorage::MCMCStorage(const MCMCParams& params)
    : n_records_((params.iter - params.n_burnin) / params.n_thin)
{
}

void MCMCStorage::initialize(const Bayes& model)
{
    mu_store_.set_size(n_records_);

    if (model.fixed().exist)
    {
        fixed_store_.set_size(model.fixed().coeff.n_elem, n_records_);
    }

    random_store_.resize(model.random().size());
    random_sigma_store_.resize(model.random().size());
    for (size_t i = 0; i < model.random().size(); ++i)
    {
        random_store_[i].set_size(model.random()[i].coeff.n_elem, n_records_);
        random_sigma_store_[i].set_size(
            model.random()[i].sigma.n_elem, n_records_);
    }

    genetic_store_.resize(model.genetic().size());
    genetic_sigma_store_.resize(model.genetic().size());
    for (size_t i = 0; i < model.genetic().size(); ++i)
    {
        genetic_store_[i].set_size(model.genetic()[i].coeff.n_elem, n_records_);
        genetic_sigma_store_[i].set_size(
            model.genetic()[i].sigma.n_elem, n_records_);
    }

    residual_store_.set_size(n_records_);
}

void MCMCStorage::store(const Bayes& model, size_t record_idx)
{
    if (record_idx >= n_records_)
    {
        return;
    }

    mu_store_[record_idx] = model.mu().value;

    if (model.fixed().exist)
    {
        fixed_store_.col(record_idx) = model.fixed().coeff;
    }

    for (size_t i = 0; i < model.random().size(); ++i)
    {
        random_store_[i].col(record_idx) = model.random()[i].coeff;
        random_sigma_store_[i].col(record_idx) = model.random()[i].sigma;
    }

    for (size_t i = 0; i < model.genetic().size(); ++i)
    {
        genetic_store_[i].col(record_idx) = model.genetic()[i].coeff;
        genetic_sigma_store_[i].col(record_idx) = model.genetic()[i].sigma;
    }

    residual_store_[record_idx] = model.residual().value;
}

}  // namespace gelex
