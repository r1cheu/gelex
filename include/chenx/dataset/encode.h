#pragma once
#include <armadillo>

namespace chenx
{
using namespace arma;
template <typename T>
Mat<T> hybird_value(
    const arma::Mat<T>& genotype,
    const arma::Col<T>& phenotype);
template <typename T>
void hybird(Mat<T>& genotype, const Mat<T>& hybird_value);
} // namespace chenx

#include "encode_impl.h"
