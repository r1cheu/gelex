// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

/**
 * @file diagnostics.h
 * @brief Diagnostics for MCMC process. refer to
 * https://github.com/pyro-ppl/numpyro/blob/master/numpyro/diagnostics.py
 */

#ifndef GELEX_BAYES_STATS_DIAGNOSTICS_H_
#define GELEX_BAYES_STATS_DIAGNOSTICS_H_
#include <Eigen/Core>
#include <utility>
#include <vector>

namespace gelex
{

// Each MatrixXd in Chains has shape (n_params, n_draws); vector length =
// n_chains.
using Chains = std::vector<Eigen::MatrixXd>;

/**
 * @brief find the smallest number >= N such that only divisor are 2, 3, 5.
 * Works just like scipy.fftpack.next_fast_len.
 * @param target N
 * @return the smallest number >= target such that only divisors are 2, 3, 5.
 */
Eigen::Index fft_next_fast_len(Eigen::Index target);

/**
 * @brief Computes R-hat over chains of samples. The samples are stored as a
 * vector of matrices where each matrix is (n_params, n_draws) and the vector
 * length is n_chains. It's required that n_chains >= 2 and n_draws >= 2.
 * Parameters that never move within a chain get NaN.
 *
 * @param samples MCMC samples
 * @return R-hat statistic for each parameter
 */
Eigen::VectorXd gelman_rubin(const Chains& samples);

/**
 * @brief Computes split R-hat over chains of samples. The samples are stored as
 * a vector of matrices where each matrix is (n_params, n_draws) and the vector
 * length is n_chains. It's required that n_draws >= 4
 *
 * @param samples
 * @return split R-hat statistic for each parameter
 */
Eigen::VectorXd split_gelman_rubin(const Chains& samples);

/**
 * @brief Compute the autocorrelation the samples at dimension n_draws
 *
 * @param x MCMC samples
 * @param bias whether to use a biased estimator
 * @return the autocorrelation of the samples
 */
Chains autocorrelation(const Chains& x, bool bias = true);

/**
 * @brief Computes the autocovariance of the samples at dimension n_draws.
 *
 * @param x MCMC Samples
 * @param bias whether to use a biased estimator
 */
Chains autocovariance(const Chains& x, bool bias = true);

/**
 * @brief Compute the effective sample size of the samples at dimension n_draws.
 * Parameters that never move get NaN.
 *
 * @param x MCMC samples, stored as a vector of matrices where each matrix is
 * (n_params, n_draws) and the vector length is n_chains
 * @param bias whether to use a biased estimator
 */
Eigen::VectorXd effect_sample_size(const Chains& x, bool bias = true);

/**
 * @brief Monte Carlo standard error of the posterior mean: the sample
 * standard deviation over every draw divided by sqrt(effective sample size).
 */
Eigen::VectorXd monte_carlo_standard_error(const Chains& x, bool bias = true);

/**
 * @brief Computes "highest posterior density interval" (HPDI) which is the
 * narrowest interval with probability mass `prob`.
 *
 * @param samples MCMC samples as a vector; left untouched
 * @param prob probability mass of the interval, in (0, 1]
 */
std::pair<double, double> hpdi(
    const Eigen::Ref<const Eigen::VectorXd>& samples,
    double prob);

auto hpdi(const Chains& chains, double prob)
    -> std::pair<Eigen::MatrixXd, Eigen::VectorXd>;
}  // namespace gelex

#endif  // GELEX_BAYES_STATS_DIAGNOSTICS_H_
