#pragma once
#include <cstddef>
#include <random>

#include "gelex/estimator/mcmc_logger.h"
#include "gelex/estimator/mcmc_params.h"
#include "gelex/estimator/mcmc_result.h"
#include "gelex/estimator/mcmc_samples.h"
#include "gelex/model/bayes.h"
#include "gelex/model/bayes_effects.h"

namespace gelex
{
class MCMC
{
   public:
    explicit MCMC(MCMCParams params);

    void run(const BayesModel& model, size_t seed = 42);
    const MCMCResult& result() const { return result_; }

   private:
    void run_one_chain(
        const BayesModel& model,
        size_t chain,
        size_t seed,
        size_t& iter);
    static void
    sample_mu(Mu& mu, arma::dvec& y_adj, double sigma_e, std::mt19937_64& rng);
    static void sample_fixed_effect(
        const FixedEffectDesign& design,
        FixedEffectState& state,
        double* y_adj,
        double sigma_e,
        std::mt19937_64& rng);
    static void sample_random_effect(
        const RandomEffectDesign& design,
        RandomEffectState& state,
        double* y_adj,
        double sigma_e,
        std::mt19937_64& rng);
    static void sample_genetic_effect(
        const GeneticEffectDesign& design,
        GeneticEffectState& state,
        double* y_adj,
        arma::uvec& snp_tracker,
        double sigma_e,
        std::mt19937_64& rng);

    MCMCLogger logger_;
    std::unique_ptr<MCMCSamples> samples_;
    MCMCParams params_;
    MCMCResult result_;
};
}  // namespace gelex
