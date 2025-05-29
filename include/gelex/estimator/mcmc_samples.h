#pragma once
#include <armadillo>
#include <vector>

#include "gelex/model/bayes.h"

namespace gelex
{

struct MCMCParams;  // Forward declaration

struct SampleGroup
{
    std::vector<arma::dmat> coeffs;
    std::vector<arma::dmat> sigmas;
};

class MCMCSamples
{
   public:
    explicit MCMCSamples(const MCMCParams& params, const BayesStatus& status);
    void store(const BayesStatus& status, size_t record_idx);

    const arma::dvec& mu() const { return mu_; }
    const arma::dmat& fixed() const { return fixed_; }
    const SampleGroup& random() const { return random_; }
    const SampleGroup& genetic() const { return genetic_; }
    const arma::dvec& residual() const { return residual_; }

    const arma::dvec& h2() const { return h2_; }

   private:
    size_t n_records_;

    arma::dvec mu_;
    arma::dmat fixed_;
    SampleGroup random_;
    SampleGroup genetic_;
    arma::dvec residual_;
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
        size_t record_idx)
    {
        for (size_t i = 0; i < status.size(); ++i)
        {
            group.coeffs[i].col(record_idx) = status[i].coeff;
            group.sigmas[i].col(record_idx) = status[i].sigma;
        }
    }
};

}  // namespace gelex
