#pragma once
#include <armadillo>
#include <cstddef>

namespace gelex
{
size_t fft_next_fast_len(size_t target);
arma::dmat gelman_rubin(const arma::dcube& samples);
arma::dmat split_gelman_rubin(const arma::dcube& samples);
arma::dcube autocorrelation(const arma::dcube& x, bool bias = true);
arma::dcube autocovariance(const arma::dcube& x, bool bias = true);
arma::dvec effect_sample_size(const arma::dcube& x, bool bias = true);
arma::dmat hpdi(const arma::dcube& samples, double prob = 0.90);
}  // namespace gelex
