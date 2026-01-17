/**
 * @file diagnostics.h
 * @brief Diagnostics for MCMC process. refer to
 * https://github.com/pyro-ppl/numpyro/blob/master/numpyro/diagnostics.py
 */

#ifndef GELEX_ESTIMATOR_BAYES_DIAGNOSTICS_H_
#define GELEX_ESTIMATOR_BAYES_DIAGNOSTICS_H_
#include <Eigen/Core>

#include "gelex/types/mcmc_samples.h"

namespace gelex
{
// the samples shape is a vector of MatrixXd. each MatrixXd is (n_params,
// n_draws). and the vector length is n_chains.

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
 * length is n_chains. It's required that n_chains >= 2 and n_draws >= 2
 *
 * @param samples MCMC samples
 * @return R-hat statistic for each parameter, shape (n_params, 1)
 */
Eigen::MatrixXd gelman_rubin(const Samples& samples);

/**
 * @brief Computes split R-hat over chains of samples. The samples are stored as
 * a vector of matrices where each matrix is (n_params, n_draws) and the vector
 * length is n_chains. It's required that n_draws >= 4
 *
 * @param samples
 * @return split R-hat statistic for each parameter, shape (n_params, 1)
 */
Eigen::MatrixXd split_gelman_rubin(const Samples& samples);

/**
 * @brief Compute the autocorrelation the samples at dimension n_draws
 *
 * @param x MCMC samples
 * @param bias whether to use a biased estimator
 * @return the autocorrelation of the samples
 */
Samples autocorrelation(const Samples& x, bool bias = true);

/**
 * @brief Computes the autocovariance of the samples at dimension n_draws.
 *
 * @param x MCMC Samples
 * @param bias whether to use a biased estimator
 */
Samples autocovariance(const Samples& x, bool bias = true);

/**
 * @brief Compute the effective sample size of the samples at dimension n_draws
 *
 * @param x MCMC samples, stored as a vector of matrices where each matrix is
 * (n_params, n_draws) and the vector length is n_chains
 * @param bias whether to use a biased estimator
 */
Eigen::VectorXd effect_sample_size(const Samples& x, bool bias = true);

/**
 * @brief Computes "highest posterior density interval" (HPDI) which is the
 narrowest interval with probability mass prob`` at dimension n_draws
 *
 * @param samples MCMC samples as a vector
 * @param prob quantiles of `samples` at `(1 - prob) / 2` and `(1 + prob) / 2`.
 */
std::pair<double, double> hpdi(
    Eigen::Ref<Eigen::VectorXd> samples,
    double prob);
}  // namespace gelex

#endif  // GELEX_ESTIMATOR_BAYES_DIAGNOSTICS_H_
