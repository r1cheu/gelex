#pragma once
#include <atomic>
#include <memory>
#include <mutex>
#include <random>

#include <Eigen/Core>

#include "../src/estimator/bayes/indicator.h"
#include "../src/logger/bayes_logger.h"
#include "gelex/estimator/bayes/params.h"
#include "gelex/estimator/bayes/result.h"
#include "gelex/estimator/bayes/samples.h"
#include "gelex/model/bayes/model.h"

namespace gelex
{
class MCMC
{
   public:
    explicit MCMC(MCMCParams params);

    const MCMCResult& run(const BayesModel& model, Eigen::Index seed = 42);
    const MCMCResult& result() const { return *result_; }

   private:
    void run_one_chain(
        const BayesModel& model,
        MCMCSamples& samples,
        Eigen::Index chain,
        Eigen::Index seed,
        std::atomic_size_t& iter,
        detail::Indicator& indicator);
    static void sample_fixed_effect(
        const bayes::FixedEffect& effect,
        bayes::FixedStatus& status,
        Eigen::Ref<Eigen::VectorXd> y_adj,
        double sigma_e,
        std::mt19937_64& rng);
    static void sample_random_effect(
        const bayes::RandomEffect& effect,
        bayes::RandomStatus& status,
        Eigen::Ref<Eigen::VectorXd> y_adj,
        double sigma_e,
        std::mt19937_64& rng);
    static void sample_additive_effect(
        const bayes::AdditiveEffect& effect,
        bayes::AdditiveStatus& status,
        Eigen::Ref<Eigen::VectorXd> y_adj,
        double sigma_e,
        std::mt19937_64& rng,
        const std::unique_ptr<GeneticTrait>& trait);

    detail::MCMCLogger logger_;

    std::mutex samples_mutex_;

    MCMCParams params_;
    std::unique_ptr<MCMCResult> result_;
};
}  // namespace gelex
