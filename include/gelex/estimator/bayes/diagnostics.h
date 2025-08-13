/**
 * @file diagnostics.h
 * @brief Diagnostics for MCMC process. refer to
 * https://github.com/pyro-ppl/numpyro/blob/master/numpyro/diagnostics.py
 */

#pragma once
#include <armadillo>
#include <cstddef>

namespace gelex
{
// the samples shape is (n_params, n_draws, n_chains) here.

/**
 * @brief find the smallest number >= N such that only divisor are 2, 3, 5.
 * Works just like scipy.fftpack.next_fast_len.
 * @param target N
 * @return the smallest number >= target such that only divisors are 2, 3, 5.
 */
size_t fft_next_fast_len(size_t target);

/**
 * @brief Computes R-hat over chains of samples, the shape of samples is
 * (n_params, n_draws, n_chains), It's required that n_chains >= 2 and n_draws
 * >= 2
 *
 * @param samples MCMC samples
 * @return R-hat statistic for each parameter, shape (n_params, 1)
 */
arma::dmat gelman_rubin(const arma::dcube& samples);

/**
 * @brief Computes split R-hat over chains of samples, the shape of samples is
 * (n_params, n_draws, n_chains), It's required that n_draws >= 4
 *
 * @param samples
 * @return split R-hat statistic for each parameter, shape (n_params, 1)
 */
arma::dmat split_gelman_rubin(const arma::dcube& samples);

/**
 * @brief Compute the autocorrelation the samples at dimension n_draws
 *
 * @param x MCMC samples
 * @param bias whether to use a biased estimator
 * @return the autocorrelation of the samples
 */
arma::dcube autocorrelation(const arma::dcube& x, bool bias = true);

/**
 * @brief Computes the autocovariance of the samples at dimension n_draws.
 *
 * @param x MCMC Samples
 * @param bias whether to use a biased estimator
 */
arma::dcube autocovariance(const arma::dcube& x, bool bias = true);

/**
 * @brief Compute the effective sample size of the samples at dimension n_draws
 *
 * @param x MCMC samples, shape (n_params, n_draws, n_chains)
 * @param bias whether to use a biased estimator
 */
arma::dvec effect_sample_size(const arma::dcube& x, bool bias = true);

/**
 * @brief Computes "highest posterior density interval" (HPDI) which is the
 narrowest interval with probability mass prob`` at dimension n_draws
 *
 * @param samples MCMC
 * @param prob quantiles of `samples` at `(1 - prob) / 2` and `(1 + prob) / 2`.
 */
std::pair<double, double> hpdi(arma::dvec& samples, double prob);
}  // namespace gelex
