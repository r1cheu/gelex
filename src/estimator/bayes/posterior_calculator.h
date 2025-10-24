#pragma once

#include <Eigen/Core>

#include "gelex/estimator/bayes/result.h"
#include "gelex/estimator/bayes/samples.h"

namespace gelex::detail
{

class EigenThreadGuard
{
   public:
    EigenThreadGuard() : old_thread_count_(Eigen::nbThreads())
    {
        Eigen::setNbThreads(1);
    }

    ~EigenThreadGuard() { Eigen::setNbThreads(old_thread_count_); }

    // Delete copy/move operations to ensure single ownership of the state.
    EigenThreadGuard(const EigenThreadGuard&) = delete;
    EigenThreadGuard& operator=(const EigenThreadGuard&) = delete;
    EigenThreadGuard(EigenThreadGuard&&) = delete;
    EigenThreadGuard& operator=(EigenThreadGuard&&) = delete;

   private:
    int old_thread_count_;
};
/**
 * @brief Posterior calculator utilities for MCMC diagnostics.
 *
 * This namespace provides statistical computations for MCMC samples,
 * eliminating code duplication and providing a clean interface.
 */
namespace PosteriorCalculator
{

/**
 * @brief Compute full posterior summary statistics (mean, std, HPDI, ESS,
 * R-hat).
 *
 * @param samples MCMC samples organized by chain (rows=params, cols=draws).
 * @param prob Probability mass for the HPDI (e.g. 0.95).
 * @return PosteriorSummary Computed posterior summary.
 */
PosteriorSummary compute_param_summary(const Samples& samples, double prob);

/**
 * @brief Compute basic posterior summary statistics (mean, std).
 *
 * @param samples MCMC samples organized by chain (rows=params, cols=draws).
 * @return PosteriorSummary Computed posterior summary with only mean and std.
 */
PosteriorSummary compute_snp_summary(const Samples& samples);

/**
 * @brief Compute mean and standard deviation for samples.
 *
 * @param summary Destination PosteriorSummary to populate with mean and std.
 * @param samples MCMC samples organized by chain (rows=params, cols=draws).
 */
void compute_mean_std(PosteriorSummary& summary, const Samples& samples);

/**
 * @brief Compute highest posterior density interval (HPDI).
 *
 * @param summary Destination PosteriorSummary to populate with HPDI.
 * @param samples MCMC samples organized by chain (rows=params, cols=draws).
 * @param prob Probability mass for the HPDI (e.g. 0.95).
 */
void compute_hpdi(
    PosteriorSummary& summary,
    const Samples& samples,
    double prob);

/**
 * @brief Compute effective sample size (ESS).
 *
 * @param summary Destination PosteriorSummary to populate with ESS.
 * @param samples MCMC samples organized by chain (rows=params, cols=draws).
 */
void compute_ess(PosteriorSummary& summary, const Samples& samples);

/**
 * @brief Compute Gelman-Rubin diagnostic (R-hat).
 *
 * @param summary Destination PosteriorSummary to populate with R-hat.
 * @param samples MCMC samples organized by chain (rows=params, cols=draws).
 */
void compute_rhat(PosteriorSummary& summary, const Samples& samples);

/**
 * @brief Compute proportion of variance explained (PVE).
 *
 * Calculates PVE = (Var(X_i) * mean(beta_i)^2) / Var(y) for each parameter.
 * Uses mean coefficients across all MCMC samples.
 *
 * @param summary Destination PosteriorSummary to populate with PVE mean.
 * stddev is set to 0.0 as PVEse is no longer calculated.
 * @param samples MCMC samples organized by chain (beta coefficients).
 * @param variances Variance of each predictor (Var(X_i)).
 * @param phenotype_var Variance of the phenotype (Var(y)).
 */
void compute_pve(
    PosteriorSummary& summary,
    const Samples& samples,
    const Eigen::VectorXd& variances,
    double phenotype_var);

/**
 * @brief Flatten multi-chain samples into a single vector for a specific
 * parameter.
 *
 * @param samples MCMC samples organized by chain.
 * @param param_index Index of the parameter to flatten.
 * @return Eigen::VectorXd Flattened samples for the parameter.
 */
Eigen::VectorXd flatten_samples(
    const Samples& samples,
    Eigen::Index param_index);

/**
 * @brief Get the number of parameters from samples.
 *
 * @param samples MCMC samples organized by chain.
 * @return Eigen::Index Number of parameters.
 */
Eigen::Index get_n_params(const Samples& samples);

/**
 * @brief Get the number of chains from samples.
 *
 * @param samples MCMC samples organized by chain.
 * @return Eigen::Index Number of chains.
 */
Eigen::Index get_n_chains(const Samples& samples);

/**
 * @brief Get the number of draws per chain from samples.
 *
 * @param samples MCMC samples organized by chain.
 * @return Eigen::Index Number of draws per chain.
 */
Eigen::Index get_n_draws(const Samples& samples);

/**
 * @brief Compute posterior inclusion probability (PIP) from SNP tracker
 * samples.
 *
 * Calculates PIP as the proportion of samples where each SNP is included
 * (tracker != 0) across all chains and draws.
 *
 * @param tracker_samples Integer samples indicating SNP inclusion status.
 * @return Eigen::VectorXd Posterior inclusion probability for each SNP.
 */
Eigen::VectorXd compute_pip(const IntSamples& tracker_samples);

}  // namespace PosteriorCalculator

}  // namespace gelex::detail
