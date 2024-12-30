// clang-format off
#pragma once
#include "chenx/kernel/gaussian_kernel.h"
#include "chenx/kernel/gram.h"

#include <armadillo>
// clang-format on

namespace chenx
{
using arma::dmat;

dmat AdditiveGrm(dmat& genotype);
dmat DomainanceGrm(dmat& genotype);
dmat ComputeGRM(dmat& genotype);

template <typename eT>
Mat<eT> Add_rbf_kernel(Mat<eT>& genotype, double bandwidth)
{
    return chenx::NaiveKernelRule<chenx::GaussianKernel, Mat<double>>::
        ApplyKernelMatrix(genotype, bandwidth);
}

template <typename eT>
Mat<eT> Dom_rbf_kernel(Mat<eT>& genotype, double bandwidth)
{
    genotype.replace(eT(2), eT(0));
    return chenx::NaiveKernelRule<chenx::GaussianKernel, Mat<double>>::
        ApplyKernelMatrix(genotype, bandwidth);
}
}  // namespace chenx
