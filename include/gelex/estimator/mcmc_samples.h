#pragma once
#include <armadillo>
#include <vector>

#include "gelex/model/bayes.h"

namespace gelex
{

struct MCMCParams;  // Forward declaration

struct SampleGroup
{
    std::vector<arma::dcube> coeffs;
    std::vector<arma::dcube> sigmas;
};

class MCMCSamples
{
   public:
    explicit MCMCSamples(const MCMCParams& params, const BayesModel& model);
    void store(const BayesStatus& status, size_t record_idx, size_t chain_idx);

    const auto& mu() const { return mu_; }
    const auto& fixed() const { return fixed_; }
    const auto& random() const { return random_; }
    const auto& genetic() const { return genetic_; }
    const auto& residual() const { return residual_; }

   private:
    size_t n_records_;
    size_t n_chains_;

    // we store mu and residule in cube(1, n_records, n_chain)
    arma::dcube mu_;
    arma::dcube fixed_;
    SampleGroup random_;
    SampleGroup genetic_;
    arma::dcube residual_;

    void init_group(
        SampleGroup& group,
        const RandomEffectDesignManager& effects) const
    {
        auto n_effects = effects.size();
        group.coeffs.resize(n_effects);
        group.sigmas.resize(n_effects);
        for (size_t i = 0; i < n_effects; ++i)
        {
            group.coeffs[i].set_size(
                effects[i].design_mat.n_cols, n_records_, n_chains_);
            group.sigmas[i].set_size(1, n_records_, n_chains_);
        }
    }

    void init_group(
        SampleGroup& group,
        const GeneticEffectDesignManager& effects) const
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

    template <typename StatusVector>
    void store_group(
        SampleGroup& group,
        const StatusVector& status,
        size_t record_idx,
        size_t chain_idx)
    {
        for (size_t i = 0; i < status.size(); ++i)
        {
            group.coeffs[i].slice(chain_idx).col(record_idx) = status[i].coeff;
            group.sigmas[i].slice(chain_idx).col(record_idx) = status[i].sigma;
        }
    }
};

}  // namespace gelex
