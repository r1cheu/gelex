#include "gelex/estimator/mcmc_diagnostics.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace gelex
{
namespace diagnostics
{

namespace
{
/**
 * @brief compute within-chain variance and variance estimator, input has shape
 * (n_params, n_draws, n_chains)
 * @param x
 * @return
 */
arma::dmat compute_chain_variance_stats(const arma::dcube& x)
{
    const size_t n_chains = x.n_slices;
    const auto n_draws = static_cast<double>(x.n_cols);
    arma::dmat chain_vars(n_chains, x.n_rows);
    arma::dmat chain_means(n_chains, x.n_rows);

    for (size_t c = 0; c < n_chains; ++c)
    {
        chain_means.row(c) = arma::mean(x.slice(c), 1).t();
        chain_vars.row(c) = arma::var(x.slice(c), 1, 1).t();
    }

    arma::dmat var_within = arma::mean(chain_vars, 0);
    arma::dmat var_estimator = var_within * (n_draws - 1) / n_draws;

    if (n_chains > 1)
    {
        arma::dmat var_between = arma::var(chain_means, 1, 0);
        var_estimator += var_between;
    }
    else
    {
        var_within = var_estimator;
    }

    return arma::join_vert(var_within, var_estimator);
}

int fft_next_fast_len(int target)
{
    if (target <= 2)
    {
        return target;
    }
    while (true)
    {
        int m = target;
        while (m % 2 == 0)
        {
            m /= 2;
        }
        while (m % 3 == 0)
        {
            m /= 3;
        }
        while (m % 5 == 0)
        {
            m /= 5;
        }
        if (m == 1)
        {
            return target;
        }
        ++target;
    }
}
}  // namespace

arma::dmat gelman_rubin(const arma::dcube& samples)
{
    arma::dmat stats = compute_chain_variance_stats(samples);
    arma::dmat rhat = arma::sqrt(stats.row(1) / stats.row(0));
    return rhat;
}

arma::dmat split_gelman_rubin(const arma::dcube& samples)
{
    const size_t N_half = samples.n_cols / 2;
    arma::dcube new_input(samples.n_rows, N_half * 2, samples.n_slices);

    for (int c = 0; c < samples.n_slices; ++c)
    {
        new_input.slice(c).cols(0, N_half - 1)
            = samples.slice(c).cols(0, N_half - 1);
        new_input.slice(c).cols(N_half, 2 * N_half - 1) = samples.slice(c).cols(
            samples.n_cols - N_half, samples.n_cols - 1);
    }

    return gelman_rubin(new_input);
}

arma::dmat autocorrelation(const arma::dmat& chain)
{
    const size_t N = chain.n_cols;
    const size_t M = fft_next_fast_len(N);
    const size_t M2 = 2 * M;

    arma::cx_mat freqvec(M2 / 2 + 1, chain.n_rows);
    arma::dmat autocorr(chain.n_rows, N);

    for (int r = 0; r < chain.n_rows; ++r)
    {
        arma::vec centered
            = arma::vectorise(chain.row(r) - arma::mean(chain.row(r)));
        arma::vec padded = arma::zeros(M2);
        padded.head(N) = centered;

        arma::cx_vec freq = arma::fft(padded);
        arma::cx_vec freq_gram = freq % arma::conj(freq);
        arma::vec autocorr_full = arma::real(arma::ifft(freq_gram));

        autocorr.row(r) = autocorr_full.head(N).t();
        autocorr.row(r) /= autocorr(r, 0);
    }

    return autocorr;
}

arma::dmat autocovariance(const arma::dmat& chain)
{
    arma::dmat autocorr = autocorrelation(chain);
    arma::vec variances = arma::var(chain, 1, 1);
    return autocorr.each_col() % variances;
}

arma::dvec effective_sample_size(const arma::dcube& samples)
{
    arma::dmat gamma_k_c(samples.n_rows, samples.n_cols);
    arma::dvec ess(samples.n_rows);

    for (int r = 0; r < samples.n_rows; ++r)
    {
        arma::dmat chain(samples.n_cols, samples.n_slices);
        for (int c = 0; c < samples.n_slices; ++c)
        {
            chain.col(c) = samples.slice(c).row(r).t();
        }
        gamma_k_c.row(r) = arma::mean(autocovariance(chain), 0);
    }

    arma::dmat stats = compute_chain_variance_stats(samples);
    arma::dmat rho_k = 1.0 - (stats.row(0) - gamma_k_c) / stats.row(1);
    rho_k.col(0).fill(1.0);

    for (int r = 0; r < rho_k.n_rows; ++r)
    {
        arma::vec Rho_k = rho_k.row(r).t();
        double tau = -1.0;

        for (int k = 0; k < Rho_k.n_elem - 1; k += 2)
        {
            double rho_sum = Rho_k(k) + Rho_k(k + 1);
            if (rho_sum < 0)
                break;
            tau += 2 * rho_sum;
        }

        ess(r) = samples.n_cols * samples.n_slices / tau;
    }

    return ess;
}

arma::dmat hpdi(const arma::dmat& samples, double prob)
{
    const int n = samples.n_rows;
    const int length = static_cast<int>(prob * n);
    arma::dmat result(2, samples.n_cols);

    for (int c = 0; c < samples.n_cols; ++c)
    {
        arma::vec sorted = arma::sort(samples.col(c));
        double min_width = std::numeric_limits<double>::max();
        int best_start = 0;

        for (int i = 0; i <= n - length; ++i)
        {
            double width = sorted(i + length - 1) - sorted(i);
            if (width < min_width)
            {
                min_width = width;
                best_start = i;
            }
        }

        result(0, c) = sorted(best_start);
        result(1, c) = sorted(best_start + length - 1);
    }

    return result;
}

}  // namespace diagnostics
}  // namespace gelex
