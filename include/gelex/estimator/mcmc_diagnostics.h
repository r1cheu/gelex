#pragma once
#include <armadillo>

namespace gelex
{

namespace diagnostics
{
arma::dmat gelman_rubin(const arma::dcube& samples);
arma::dmat split_gelman_rubin(const arma::dcube& samples);
arma::dvec effective_sample_size(const arma::dcube& samples);
arma::dmat hpdi(const arma::dmat& samples, double prob = 0.90);
arma::dmat autocorrelation(const arma::dmat& chain);
arma::dmat autocovariance(const arma::dmat& chain);
}  // namespace diagnostics
}  // namespace gelex
