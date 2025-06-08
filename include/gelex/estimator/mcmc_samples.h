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
    explicit MCMCSamples(const MCMCParams& params, const BayesStatus& status, size_t n_chains);
    void store(const BayesStatus& status, size_t record_idx, size_t chain_idx);

    const arma::dmat& mu() const { return mu_; }
    const arma::dcube& fixed() const { return fixed_; }
    const SampleGroup& random() const { return random_; }
    const SampleGroup& genetic() const { return genetic_; }
    const arma::dmat& residual() const { return residual_; }

    const arma::dvec& h2() const { return h2_; }

   private:
    size_t n_records_;
    size_t n_chains_;

    arma::dmat mu_;
    arma::dcube fixed_;
    SampleGroup random_;
    SampleGroup genetic_;
    arma::dmat residual_;
    arma::dvec h2_;

    template <typename StatusVector>
    void init_group(SampleGroup& group, const StatusVector& status)
    {
        group.coeffs.resize(status.size());
        group.sigmas.resize(status.size());
        for (size_t i = 0; i < status.size(); ++i)
        {
            group.coeffs[i].set_size(status[i].coeff.n_elem, n_records_);
            group.sigmas[i].set_size(status[i].sigma.n_elem, n_records_);
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
