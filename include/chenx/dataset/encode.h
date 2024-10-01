#pragma once
#include <armadillo>

namespace chenx {
namespace dataset {
template <typename T>
void dominance(arma::Mat<T>& genotype);
template <typename T>
arma::Mat<T> hybird_value(const arma::Mat<T>& genotype, const arma::Col<T>& phenotype);
template <typename T>
void hybird(arma::Mat<T>& genotype, const arma::Mat<T>& hybird_value);
}  // namespace dataset
}  // namespace chenx

#include "encode_impl.h"
