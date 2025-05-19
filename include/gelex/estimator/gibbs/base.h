#pragma once

#include <armadillo>
#include <cmath>

#include "gelex/dist.h"

namespace gelex
{
#ifdef ARMA_USE_LAPACK
#if !defined(ARMA_BLAS_CAPITALS)
#define arma_daxpy daxpy
#else
#define arma_daxpy DAXPY
#endif
extern "C" void arma_fortran(arma_daxpy)(
    const int* n,
    const double* da,
    const double* dx,
    const int* incx,
    double* dy,
    const int* incy);
#endif

using dvec = arma::vec;

inline void daxpy_auto(arma::dvec& y, const arma::dvec& x, double alpha)
{
    const int n = x.n_elem;
    static const int inc = 1;

    if (n <= 10000)
    {
        y += alpha * x;
    }
    else
    {
        daxpy_(&n, &alpha, x.memptr(), &inc, y.memptr(), &inc);
    }
}

inline double
compute_rhs(const dvec& col_i, const dvec& y_adj, double old_i, double col_norm)
{
    return arma::dot(col_i, y_adj) + (col_norm * old_i);
}

inline void sample_effect(
    Normal& normal,
    dvec& coeff,
    dvec& y_adj,
    const arma::dmat& design_mat,
    const dvec& cols_norm2,
    double sigma_e,
    double sigma)
{
    const dvec inv_scaler = 1.0 / (cols_norm2 + sigma_e / sigma);

    for (size_t i = 0; i < coeff.n_elem; ++i)
    {
        const double old_i = coeff.at(i);
        const dvec& col_i = design_mat.unsafe_col(i);
        const double inv_scaler_i = inv_scaler.at(i);
        double rhs = compute_rhs(col_i, y_adj, old_i, cols_norm2.at(i));
        double new_i = normal(rhs * inv_scaler_i, sqrt(sigma_e * inv_scaler_i));
        coeff.at(i) = new_i;
        daxpy_auto(y_adj, col_i, old_i - new_i);
    }
}

}  // namespace gelex
