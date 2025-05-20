#pragma once
#include <memory>
#include <random>

#include "gelex/estimator/mcmc_logger.h"
#include "gelex/estimator/mcmc_params.h"
#include "gelex/estimator/mcmc_result.h"
#include "gelex/estimator/mcmc_storage.h"
#include "gelex/model/bayes.h"
#include "gelex/model/effects/bayes_effects.h"

namespace gelex
{
class MCMC
{
   public:
    explicit MCMC(MCMCParams params);

    const MCMCResult& run(Bayes& model, size_t log_freq);
    const MCMCResult& result() const { return result_; }
    const MCMCStorage& storage() const { return *storage_; }

   private:
    void sample_mu(Bayes& model);
    void sample_fixed_effect(bayes::FixedEffect& effect, double sigma_e);
    void sample_random_effect(bayes::RandomEffect& effect, double sigma_e);
    void sample_genetic_effect(bayes::GeneticEffect& effect, double sigma_e);

    MCMCLogger logger_;

    std::mt19937_64 rng_;
    arma::dvec y_adj_;
    arma::uvec snp_tracker_;
    MCMCParams params_;
    std::unique_ptr<MCMCStorage> storage_;
    MCMCResult result_;
};
}  // namespace gelex
