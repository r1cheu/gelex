#pragma once

#include <armadillo>
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

inline void daxpy_ptr(int n, double alpha, const double* x, double* y)
{
    static const int inc = 1;
    daxpy_(&n, &alpha, x, &inc, y, &inc);
}

inline double ddot_ptr(int n, const double* x, const double* y)
{
    static const int inc = 1;
    return arma::ddot_(&n, x, &inc, y, &inc);
}

}  // namespace gelex
