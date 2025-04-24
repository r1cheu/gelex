#pragma once
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <armadillo>

namespace gelex
{
using arma::dmat;
using arma::dvec;
using arma::uvec;

struct sigma_prior
{
    double nu;
    double s2;
};

struct Priors
{
    sigma_prior sigma_a{-2, 0};
    sigma_prior sigma_r{2, 0};
    sigma_prior sigma_e{-2, 0};
    dvec pi;
};
}  // namespace gelex
