#pragma once
#include <armadillo>
#include "chenx/kernel/gaussian_kernel.h"
#include "chenx/kernel/gram.h"

namespace chenx
{
using namespace arma;
template <typename eT>
Mat<eT> Amat(Mat<eT>& genotype)
{
    Col<eT> pA = mean(genotype, 1) / 2;
    genotype.each_col() -= 2 * pA;
    Mat<eT> grm = genotype.t() * genotype;
    return grm / trace(grm) * grm.n_cols;
}

template <typename eT>
Mat<eT> Dmat(Mat<eT>& genotype)
{
    Col<eT> pA = mean(genotype, 1) / 2;
    Col<eT> pa = 1 - pA;
    genotype.replace(eT(2), eT(0));
    genotype.each_col() -= 2 * (pA % pa);
    Mat<eT> grm = genotype.t() * genotype;
    return grm / trace(grm) * grm.n_cols;
}

template <typename eT>
Mat<eT> Add_rbf_kernel(Mat<eT>& genotype, double bandwidth)
{
    return chenx::NaiveKernelRule<chenx::GaussianKernel, Mat<double>>::ApplyKernelMatrix(genotype, bandwidth);
}

template <typename eT>
Mat<eT> Dom_rbf_kernel(Mat<eT>& genotype, double bandwidth)
{
    genotype.replace(eT(2), eT(0));
    return chenx::NaiveKernelRule<chenx::GaussianKernel, Mat<double>>::ApplyKernelMatrix(genotype, bandwidth);
}
}  // namespace chenx
