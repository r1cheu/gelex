#pragma once
#include <armadillo>

namespace chenx {
namespace dataset {

template <typename T>
void mean_impute(arma::Mat<T>& genotype);

template <typename T>
void median_impute(arma::Mat<T>& genotype);
}  // namespace dataset
}  // namespace chenx

#include "impute_impl.h"
