#pragma once
#include <armadillo>
#include <vector>

#include "gelex/model/bayes.h"

namespace gelex
{

struct MCMCParams;  // Forward declaration

class MCMCStorage
{
   public:
    explicit MCMCStorage(const MCMCParams& params);

    void initialize(const Bayes& model);
    void store(const Bayes& model, size_t record_idx);

    const arma::dvec& mu_samples() const { return mu_store_; }
    const arma::dmat& fixed_samples() const { return fixed_store_; }
    const std::vector<arma::dmat>& random_samples() const
    {
        return random_store_;
    }
    const std::vector<arma::dmat>& genetic_samples() const
    {
        return genetic_store_;
    }
    const arma::dvec& residual_samples() const { return residual_store_; }
    const std::vector<arma::dmat>& random_sigma_samples() const
    {
        return random_sigma_store_;
    }
    const std::vector<arma::dmat>& genetic_sigma_samples() const
    {
        return genetic_sigma_store_;
    }

   private:
    size_t n_records_;

    // Effect samples
    arma::dvec mu_store_;
    arma::dmat fixed_store_;
    std::vector<arma::dmat> random_store_;
    std::vector<arma::dmat> genetic_store_;
    arma::dvec residual_store_;

    // Variance components
    std::vector<arma::dmat> random_sigma_store_;
    std::vector<arma::dmat> genetic_sigma_store_;
};

}  // namespace gelex
