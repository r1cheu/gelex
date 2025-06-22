#include "gelex/estimator/bayes/diagnostics.h"

#include <cstddef>

#include <armadillo>

namespace gelex
{

using arma::dcube;
using arma::dmat;
using arma::dvec;
/**
 * @brief compute within-chain variance and variance estimator, input has shape
 * (n_params, n_draws, n_chains)
 * @param x
 * @return
 */
std::pair<dvec, dvec> compute_chain_variance_stats(const dcube& x)
{
    const size_t n_chains = x.n_slices;
    const auto n_draws = static_cast<double>(x.n_cols);

    dmat chain_vars(x.n_rows, n_chains);
    dmat chain_means(x.n_rows, n_chains);

    for (size_t c = 0; c < n_chains; ++c)
    {
        chain_means.col(c) = arma::mean(x.slice(c), 1);
        chain_vars.col(c) = arma::var(x.slice(c), 0, 1);
    }

    dvec var_within = arma::mean(chain_vars, 1);
    dvec var_estimator = var_within * (n_draws - 1) / n_draws;

    if (n_chains > 1)
    {
        dmat var_between = arma::var(chain_means, 0, 1);
        var_estimator += var_between;
    }
    else
    {
        var_within = var_estimator;
    }

    return {var_within, var_estimator};
}

size_t fft_next_fast_len(size_t target)
{
    if (target <= 2)
    {
        return target;
    }
    while (true)
    {
        size_t m = target;
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

dmat gelman_rubin(const dcube& samples)
{
    auto [var_within, var_estimator] = compute_chain_variance_stats(samples);
    dmat rhat = arma::sqrt(var_estimator / var_within);
    return rhat;
}

dmat split_gelman_rubin(const dcube& samples)
{
    const size_t N_half = samples.n_cols / 2;
    dcube new_input = arma::join_slices(
        samples.cols(0, N_half - 1), samples.cols(N_half, samples.n_cols - 1));

    return gelman_rubin(new_input);
}

arma::mat autocorrelation(const dmat& x, bool bias)
{
    dmat signal = x.t();
    const size_t N = signal.n_rows;
    const size_t M = fft_next_fast_len(N);
    const size_t M2 = 2 * M;

    signal.each_col([](dvec& col) { col -= arma::mean(col); });

    arma::cx_mat freqvec = arma::fft(signal, M2);
    arma::cx_mat freqvec_gram = freqvec % arma::conj(freqvec);

    arma::cx_mat autocorr_cx = arma::ifft(freqvec_gram, M2);
    dmat autocorr = arma::real(autocorr_cx.head_rows(N));

    if (!bias)
    {
        autocorr.each_col() /= arma::regspace(N, -1, 1);
    }

    arma::rowvec variances = autocorr.row(0);
    autocorr.each_row() /= variances;

    return autocorr.t();
}

dcube autocorrelation(const dcube& x, bool bias)
{
    dcube result(x.n_rows, x.n_cols, x.n_slices);
    for (size_t s = 0; s < x.n_slices; ++s)
    {
        result.slice(s) = autocorrelation(x.slice(s), bias);
    }
    return result;
}

dcube autocovariance(const dcube& x, bool bias)
{
    dcube result = autocorrelation(x, bias);

    dcube x_var(x.n_rows, 1, x.n_slices);

    for (size_t s = 0; s < x.n_slices; ++s)
    {
        x_var.slice(s) = arma::var(x.slice(s), 1, 1);
    }

    for (size_t s = 0; s < result.n_slices; ++s)
    {
        result.slice(s).each_col([&x_var, s](dvec& col)
                                 { col %= x_var.slice(s); });
    }

    return result;
}

dvec effect_sample_size(const dcube& x, bool bias)
{
    const size_t n_params = x.n_rows;
    const size_t n_draws = x.n_cols;
    const size_t n_chains = x.n_slices;

    if (n_draws < 2)
    {
        throw std::invalid_argument("At least 2 draws are required");
    }

    dcube gamma_k_c = autocovariance(x, bias);
    dmat gamma_k_c_mean = arma::mean(gamma_k_c, 2);
    auto [var_within, var_estimator] = compute_chain_variance_stats(x);

    dmat var_within_boardcast
        = var_within * arma::ones<arma::rowvec>(1, n_draws);
    dmat var_estimator_boardcast
        = var_estimator * arma::ones<arma::rowvec>(1, n_draws);
    dmat rho_k(n_params, n_draws, arma::fill::ones);
    rho_k -= (var_within_boardcast - gamma_k_c_mean) / var_estimator_boardcast;
    rho_k.col(0).fill(1);

    const size_t n_pairs = n_draws / 2;
    dmat Rho_k(n_params, n_pairs);

    for (size_t j = 0; j < n_pairs; ++j)
    {
        Rho_k.col(j) = rho_k.col(2 * j) + rho_k.col((2 * j) + 1);
    }

    dmat Rho_mono = Rho_k;

    for (size_t i = 0; i < n_params; ++i)
    {
        double current_min = Rho_k(i, 0);

        for (size_t j = 1; j < n_pairs; ++j)
        {
            double val = std::max(Rho_k.at(i, j), 0.0);
            val = std::min(val, current_min);
            current_min = val;
            Rho_mono.at(i, j) = val;
        }
    }

    dvec Rho_sum = arma::sum(Rho_mono, 1);
    dvec s2 = -1.0 + 2.0 * Rho_sum;
    auto total_samples = static_cast<double>(n_chains * n_draws);
    dvec n_eff = total_samples / s2;

    return n_eff;
}

arma::rowvec hpdi(arma::dvec samples, double prob)
{
    std::sort(samples.begin(), samples.end());
    if (prob == 1)
    {
        return {samples.front(), samples.back()};
    }
    size_t mass = samples.n_elem;
    auto index_length = static_cast<size_t>(prob * static_cast<double>(mass));
    size_t tails = mass - index_length;

    dvec intervals_left = samples.head(tails);
    dvec intervals_right = samples.tail(tails);

    dvec intervals = intervals_right - intervals_left;

    size_t index_start = intervals.index_min();
    size_t index_end = index_start + index_length;

    return {samples.at(index_start), samples.at(index_end)};
}

arma::dmat hpdi(const arma::dcube& samples, double prob)
{
    if (prob <= 0.0 || prob > 1.0)
    {
        throw std::invalid_argument("Probability must be in (0, 1]");
    }
    arma::dmat result(samples.n_rows, 2);
#pragma omp parallel for default(none) shared(samples, result, prob)
    for (size_t r = 0; r < samples.n_rows; ++r)
    {
        result.row(r) = hpdi(arma::vectorise(samples.row(r)), prob);
    }
    return result;
}

}  // namespace gelex
