#pragma once
#include <cstddef>
#include <memory>
#include <mutex>
#include <random>

#include "gelex/estimator/bayes/indicator.h"
#include "gelex/estimator/bayes/logger.h"
#include "gelex/estimator/bayes/params.h"
#include "gelex/estimator/bayes/result.h"
#include "gelex/estimator/bayes/samples.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/model.h"

namespace gelex
{
class MCMC
{
   public:
    explicit MCMC(MCMCParams params);

    const MCMCResult& run(const BayesModel& model, size_t seed = 42);
    const MCMCResult& result() const { return *result_; }
    const MCMCSamples& samples() const { return *samples_; }

   private:
    void run_one_chain(
        const BayesModel& model,
        size_t chain,
        size_t seed,
        std::atomic_size_t& iter,
        Indicator& indicator);
    static void sample_fixed_effect(
        const bayes::FixedEffect& design,
        bayes::FixedEffectState& state,
        double* y_adj,
        double sigma_e,
        std::mt19937_64& rng);
    static void sample_random_effect(
        const bayes::RandomEffect& design,
        bayes::RandomEffectState& state,
        double* y_adj,
        double sigma_e,
        std::mt19937_64& rng);
    static void sample_genetic_effect(
        const bayes::GeneticEffect& design,
        bayes::GeneticEffectState& state,
        double* y_adj,
        arma::uvec& snp_tracker,
        double sigma_e,
        std::mt19937_64& rng);

    MCMCLogger logger_;

    std::unique_ptr<MCMCSamples> samples_;
    std::mutex samples_mutex_;

    MCMCParams params_;
    std::unique_ptr<MCMCResult> result_;
};
}  // namespace gelex
