#pragma once

#include <armadillo>

namespace gelex
{

namespace detail
{
std::pair<double, double> compute_chain_variance_stats(const arma::dmat&);
std::pair<double, double> compute_chain_variance_stats(const arma::dcube&);
}  // namespace detail
}  // namespace gelex
